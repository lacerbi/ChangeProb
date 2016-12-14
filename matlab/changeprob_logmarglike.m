function [lmL, modelPost, nLL, rmse, fitParams, resp_model,...
    resp_obs, p_true, p_estimate, post] = changeprob_logmarglike(model, data, task, parameters, gridSize, paramBounds, fixNoise)

%% CHANGEPROB_LOGMARGLIKE Log marginal likelihood for specified model

% Here we compute the log marginal likelihood of the specfied model. 
% To calculate the marginal likelihood, we multiply the probability of 
% the data given the model by the prior of the model, integrate over 
% the model parameters, and then take the logarithm.

% INPUT:
    % model: a character string indicating the model you want to fit, 
    % which is limited to the following:
        % 'idealBayesian'
        % 'fixed'
        % 'exponential'
        % 'RL_probability'
        % 'RL_criterion'
    % data: data struct from changing probability experiment
    % task: lets you choose which task to fit
        % 1 - overt-criterion task
        % 2 - covert-criterion task
        % 3 - mixed design (overt-criterion task on every 5th trial)
    % Vector indicating the parameters to-be-fit (1 - fit, 0 - not fit)
        % [sigma_ellipse, sigma_criterion, lapse, gamma, alpha, w, Tmax]
    % gridSize: size of parameter grid (e.g., [n x m x p])
    % paramBounds: lower and upper parameter bounds (numel(gridSize) x 2)
    % fixNoise: fix noise from measurement session (leave empty for default,
    %           which is 0 for covert, 1 for overt)

% OUTPUT:
    % lmL: a measure of model fit
    % modelPost: proportional to the p(parameters | data, model)
    % nLL: negative log likelihood that correponds to the best fitting
    % parameters
    % rmse: root mean squared error between the model probability and
    % the true probability (Note: does not apply to the RL_criterion model)
    % fitParams: best fitting model parameters computed by taking the
    % MAP of modelPost
    % resp_model: predicted criterion
    % resp_obs: observer's criterion
    % p_true: true probability values
    % p_estimate: model's estimate of probability (Note: does not apply to
    % the RL_criterion model)

% Authors:  Elyse Norton, Luigi Acerbi
% Email:    {elyse.norton,luigi.acerbi}@gmail.com
% Date:     12/14/2016
    
%tic
% Model to be fit
if nargin < 1; error('Please indicate the model you want to fit.'); end
potentialModels = {'idealBayesian', 'fixed', 'exponential', 'RL_probability', ...
    'RL_criterion'};
model = find(strcmp(model, potentialModels)); % recode model to numeric value

% Data struct or random seed for fake data generation
if nargin < 2 || isempty(data); data = 0; end
if isnumeric(data); rng(data); data = []; end

% Task (1 overt, 2 covert)
if nargin < 3 || isempty(task); task = 1; end % overt-criterion task is the default
if task ~= 1 && task ~= 2 && task ~= 3; error('TASK can only be 1 (overt), 2 (covert), or 3 (mixed).'); end
if task == 1; taskName = 'overt'; elseif task == 2; taskName = 'covert'; else taskName = 'mixed'; end

NumSamples = 5000;
MaxParams = 7;

% Parameters to be fit
if nargin < 4 || isempty(parameters)
    parameters = zeros(1,MaxParams);    
    switch task
        case 1
            if model < 3
                parameters(2) = 1;
            else
                parameters([2,5]) = 1;
            end
        case 2
            if model < 3
                parameters(1) = 1;
            else
                parameters([1,5]) = 1;
            end
        case 3
            if model < 3
                parameters(2) = 1;
            else
                parameters([2,5]) = 1;
            end
    end
    if nargin < 4
        gridSize = [];
        paramBounds = [];
    end
elseif sum(parameters) == 0
    error('You must specify at least one parameter to fit.')
end

if nargin < 7
     fixNoise = [];
end

if isempty(fixNoise)
    switch task
        case 1; parameters(1) = 0;  % Overt task: by default noise is fixed
        case 2; parameters(1) = 1;  % Covert task: by default noise is free
        case 3; parameters(1) = 0;  % Mixed task: by default noise is fixed
    end
else
    if fixNoise; parameters(1) = 0; else parameters(1) = 1; end
end

paramNames = {'sigma_ellipse', 'sigma_criterion', 'lambda', 'gamma', 'alpha', 'w', 'Tmax'};
NumParams = sum(parameters);
if NumParams > 0; fitParamNames = paramNames{logical(parameters)}; end 
I_params = find(parameters ~= 0);
% Sampling rate of parameter grid
if nargin < 5 || isempty(gridSize)
    gridSize = 100*ones(1, NumParams); % Default grid is 100 x 100
    if nargin < 5
        paramBounds = [];
    end
elseif numel(gridSize) ~= NumParams
    error('Matrix dimensions do not agree. numel(gridSize) must equal sum(parameters).');
end

% Lower and upper parameter bounds
if nargin < 6 || isempty(paramBounds)    
    paramBounds_def = [1,30; 1,30; 0,0.1; -Inf,Inf; 0,1; 0,1; 2,200];    
    paramBounds = paramBounds_def(I_params,:);
end

if NumParams ~= size(paramBounds,1)
    error('Please specify parameter bounds for each parameter to be fit.');
end

%% Print session
if fixNoise; noiseString = 'FIXED'; else noiseString = 'FREE'; end
fprintf('Fitting model %s, %s-criterion task; sensory noise is %s, %d free parameters.\n', potentialModels{model}, taskName, noiseString, NumParams);

%% Get session parameters
[NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs, score] = changeprob_getSessionParameters(data, task);
if and(task == 1, model == 5) || and(task == 3, model == 5)
    X = bsxfun(@plus, S, sigma_ellipse*randn(numel(S), NumSamples));
else
    X = [];
end

%% Observer model parameters
if NumParams > 0
    params2fit = zeros(NumParams, gridSize(1));
    for iParam = 1:NumParams
        switch I_params(iParam)
            case {1, 2}
                params2fit(iParam,:) = log(linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam))); % sigma_noise
            case 3
                params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % lambda
            case 4
                params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % gamma
            case 5
                params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % alpha
            case 6
                params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % w
            case 7
                params2fit(iParam,:) = round(linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam))); % Tmax (discrete)
        end
    end

    % Create parameter list for fitting
    for iParam = 1:NumParams; params2fit_cell{iParam} = params2fit(iParam,:); end
    params2fit_list = combvec(params2fit_cell{:})';    
    exp_idx = I_params == 1 | I_params == 2;
    params2fit_list(:,exp_idx) = exp(params2fit_list(:,exp_idx));   % Check this
    clear params2fit_cell;
end

% Default parameter values for those not fitted
inputParams = zeros(1,numel(parameters));
I_notFit = find(parameters == 0);
sigmacriterion_def = 5; % Default sigma_criterion
lapse_def = 1e-4;       % Default lapse (i.e., tiny lapse)
gamma_def = Inf;        % Default gamma (i.e., BDT)
alpha_def = 0.2;        % Default alpha
w_def     = 1;          % Default w (i.e., no bias)
Tmax_def  = 0;          % None provided, use default prior window
notFit_def = [sigma_ellipse, sigmacriterion_def, lapse_def, gamma_def, alpha_def, w_def, Tmax_def];
inputParams(I_notFit) = notFit_def(I_notFit);

inputParams

%% Separate case if there are no free parameters
if NumParams == 0
    [nLL,rmse,resp_model,p_estimate,post] = ...
        changeprob_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, task, score, model, X);
    lmL = -nLL;     % Log marginal likelihood is simply log likelihood
    modelPost = [];
    fitParams = [];    
    return;
end

%% Choose priors - start with uniformative priors for all parameters
    % These will be uniform for all priors since we took the log of the
    % sigma parameters
    
% Uniform prior for all parameters
prior = 1/prod(params2fit(:,end) - params2fit(:,1));
logprior = log(prior);

%% Compute the negative log likelihood for all parameter sets

nLL_mat = zeros([gridSize,1]);
maxii = prod(gridSize);
for ii = 1:maxii
    if rem(ii,500) == 0; fprintf('%.1f%%..', 100*ii/maxii); end
    % Parameters in params2fit_list are already transformed
    inputParams(I_params) = params2fit_list(ii,:);
    nLL_mat(ii) = changeprob_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, task, score, model, X);
end
fprintf('\n');
clear params2fit_list;  % Free up some memory

%% Compute the marginal likelihood

% Add the log likelihood to the log prior, subtract the max, and exponentiate
logUpost = bsxfun(@plus, -nLL_mat, logprior);
maxlogUpost = max(logUpost(:));
modelUpost = exp(logUpost - maxlogUpost);   % Unnormalized posterior       

% Trapezoidal integration across all parameters
temp = modelUpost;
for iParam = 1:NumParams; temp = qtrapz(temp); end
dV = prod(params2fit(:,2) - params2fit(:,1));   % volume element
Z = temp*dV;    % Normalization constant

% Normalized posterior
modelPost = modelUpost / Z;

% Log marginal likelihood (add back max logUpost from before)
lmL = log(Z) + maxlogUpost;

% Find the MAP (just for now; we should do better)
[~,linidx] = max(logUpost(:));
[idx(1),idx(2),idx(3),idx(4),idx(5)] = ind2sub(size(logUpost),linidx);
for iParam = 1:NumParams; bestFit_param(iParam) = params2fit(iParam,idx(iParam)); end

% Transform parameters (check that this is correct)
fitParams = bestFit_param;
fitParams(exp_idx) = exp(fitParams(exp_idx));

% Recompute additional outputs from best fit params
inputParams(I_params) = fitParams;
[nLL,rmse,resp_model,p_estimate,post] = ...
    changeprob_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, task, score, model, X);

%toc
end