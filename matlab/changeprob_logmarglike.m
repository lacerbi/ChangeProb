function [lmL, modelPost, nLL, rmse, fitParams, resp_model,...
    resp_obs, p_true, p_estimate, post] = changeprob_logmarglike(model, data, task, parameters, gridSize, paramBounds)

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
    % Vector indicating the parameters to-be-fit (1 - fit, 0 - not fit)
        % [sigma_ellipse, sigma_criterion, lapse, gamma, alpha, w]
    % gridSize: size of parameter grid (e.g., [n x m x p])
    % paramBounds: lower and upper parameter bounds (numel(gridSize) x 2) 

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
% Date:     11/06/2016
    
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
if task ~= 1 && task ~= 2; error('TASK can only be 1 (overt) or 2 (covert).'); end

% Parameters to be fit
if nargin < 4 || isempty(parameters)
    switch task
        case 1
            if model < 3
                parameters = [0 1 0 0 0 0];
            else
                parameters = [0 1 0 0 1 0];
            end
        case 2
            if model < 3
                parameters = [1 0 0 0 0 0];
            else
                parameters = [1 0 0 0 1 0];
            end
    end
    if nargin < 4
        gridSize = [];
        paramBounds = [];
    end
elseif sum(parameters) == 0
    error('You must specify at least one parameter to fit.')
end

paramNames = {'sigma_ellipse', 'sigma_criterion', 'lambda', 'gamma', 'alpha', 'w'};
fitParamNames = paramNames{logical(parameters)};
NumParams = sum(parameters);
I_params = find(parameters ~= 0);
% Sampling rate of parameter grid
if nargin < 5 || isempty(gridSize)
    gridSize = zeros(1, NumParams) + 100; % Default grid is 100 x 100
    if nargin < 5
        paramBounds = [];
    end
elseif numel(gridSize) ~= NumParams
    error('Matrix dimensions do not agree. numel(gridSize) must equal sum(parameters).');
end

% Lower and upper parameter bounds
if nargin < 6 || isempty(paramBounds)    
    paramBounds_def = [1,30; 1,30; 0,0.1; -Inf,Inf; 0,1; 0,1];    
    paramBounds = paramBounds_def(I_params,:);
end

if NumParams ~= size(paramBounds,1)
    error('Please specify parameter bounds for each parameter to be fit.');
end

%% Get session parameters
[NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs, score] = changeprob_getSessionParameters(data, task);
if task == 1 && model == 5
    X = bsxfun(@plus, S, sigma_ellipse*randn(numel(S), 1000));
else
    X = zeros(numel(S), 1000);
end

%% Observer model parameters
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
    end
end

% Default parameter values for those not fitted
inputParams = zeros(1,numel(parameters));
I_notFit = find(parameters == 0);
sigmacriterion_def = 5; % Default sigma_criterion
lapse_def = 0;          % Default lapse (i.e., no lapse)
gamma_def = Inf;        % Default gamma (i.e., BDT)
alpha_def = 0.2;        % Default alpha
w_def     = 1;          % Default w (i.e., no bias)
notFit_def = [sigma_ellipse, sigmacriterion_def, lapse_def, gamma_def, alpha_def, w_def];
inputParams(I_notFit) = notFit_def(I_notFit);

%% Choose priors - start with uniformative priors for all parameters
    % These will be uniform for all priors since we took the log of the
    % sigma parameters
    
% Uniform prior for all parameters    
prior = 1/prod(params2fit(:,end) - params2fit(:,1));
logprior = log(prior);

%% Compute the negative log likelihood for all parameter sets

nLL_mat = zeros(gridSize);
fitParams = zeros(1,NumParams);
maxii = prod(gridSize);
for ii = 1:maxii
    if rem(ii,500) == 0; fprintf('%.1f%%..', 100*ii/maxii); end
    [idx(1),idx(2),idx(3),idx(4),idx(5)] = ind2sub(gridSize,ii);
    for iParam = 1:NumParams; fitParams(iParam) = params2fit(iParam,idx(iParam)); end
    fitParams(1) = exp(fitParams(1));
    inputParams(I_params) = fitParams;
    nLL_mat(ii) = changeprob_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, task, score, model, X);
end
fprintf('\n');

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
for iParam = 1:NumParams; bestFit_param(iParam) = params2fit(idx(iParam)); end

% Transform parameters (check that this is correct)
fitParams = bestFit_param; fitParams(1) = exp(fitParams(1));

% Recompute additional outputs from best fit params
inputParams(I_params) = fitParams;
[nLL,rmse,resp_model,p_estimate,post] = ...
    changeprob_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, task, score, model, X);

%toc
end