function [nLL, fitParams, resp_model, resp_obs, p_true, p_estimate] = changeprob_maxlike(model,...
    data, task, parameters, paramBounds)
%% CHANGEPROB_MAXLIKE Max log likelihood for specified model

% [Documentation to be written]

if nargin < 2; data = []; end
if nargin < 3; task = []; end
if nargin < 4; parameters = []; end
if nargin < 5; paramBounds = []; end

% Model to be fit
if nargin < 1; error('Please indicate the model you want to fit.'); end
potentialModels = {'idealBayesian', 'fixed', 'exponential', 'RL_probability', ...
    'RL_criterion', 'gold'};
model = find(strcmp(model, potentialModels)); % recode model to numeric value

% Data struct or random seed for fake data generation
if isempty(data); data = 0; end
if isnumeric(data); rng(data); data = []; end

% Task (1 overt, 2 covert)
if isempty(task); task = 1; end % overt-criterion task is the default
if task ~= 1 && task ~= 2 && task ~= 3; error('TASK can only be 1 (overt), 2 (covert), or 3 (mixed).'); end
if task == 1; taskName = 'overt'; elseif task == 2; taskName = 'covert'; else taskName = 'mixed'; end

NumSamples = 5000;
MaxParams = 15;

% Parameters to be fit
if isempty(parameters)
    parameters = zeros(1,MaxParams);    
    switch task
        case 1
            if model < 3
                parameters(2) = 1;
            elseif and(model > 2, model ~=6)
                parameters([2,5]) = 1;
            else
                parameters([2,10,11]) = 1;
            end
        case 2
            if model < 3
                parameters(1) = 1;
            elseif and(model > 2, model ~=6)
                parameters([1,5]) = 1;
            else
                parameters([1,10,11]) = 1;
            end
        case 3
            if model < 3
                parameters(2) = 1;
            elseif and(model > 2, model ~=6)
                parameters([2,5]) = 1;
            else
                parameters([2,10,11]) = 1;
            end
    end
elseif sum(parameters) == 0
    error('You must specify at least one parameter to fit.')
end

if numel(parameters) < MaxParams
    for idx_padding = 1:(MaxParams - numel(parameters))
        parameters(:,end+1) = 0;
    end
end

paramNames = {'sigma_ellipse', 'sigma_criterion', 'lambda', 'gamma', 'alpha', 'w', 'Tmax', ...
    'pVec', 'beta', 'delta1', 'delta2', 'hRate', 'nu_p', 'delta_Tmax', 'delta_pVec'};
NumParams = sum(parameters);
I_params = find(parameters ~= 0);

% Lower and upper parameter bounds
if isempty(paramBounds)    
    paramBounds_def = [1,30; 1,30; 0,0.1; -Inf,Inf; 0,1; 0,1; 2,200; 0,.5; 0,10; 1.01,5; 1.01,14; 0,1; 0,5; 2,200; 0,.5];    
    paramBounds = paramBounds_def(I_params,:);
end

if NumParams ~= size(paramBounds,1)
    error('Please specify parameter bounds for each parameter to be fit.');
end

%% Get session parameters

[NumTrials, sigma_ellipseData, mu, sigma_s, C, S, p_true, resp_obs, score] = changeprob_getSessionParameters(data, task);

if and(task == 1, model == 5) || and(task == 3, model == 5)
    X = bsxfun(@plus, S, sigma_ellipseData*randn(numel(S), NumSamples));
else
    X = [];
end

% Default parameter values for those not fitted
inputParams = zeros(1,numel(parameters));
I_notFit = find(parameters == 0);
sigmacriterion_def = 5; % Default sigma_criterion
lapse_def = 1e-4;       % Default lapse (i.e., tiny lapse)
gamma_def = Inf;        % Default gamma (i.e., BDT)
alpha_def = 0.2;        % Default alpha
w_def     = 1;          % Default w (i.e., no bias)
Tmax_def  = 0;          % Use default prior window
pVec_def = 0;           % Use default probability vector
beta_def = 0;           % Use default hyperprior, [0,0]
delta_def = 2;          % Use default node distance
hRate_def = .01;        % Use default hazard rate (average rate of change)
nu_p_def = log(2);      % Use default nu_p (Beta(1,1))
delta_Tmax_def = 0;     % Use default delta_Tmax
delta_pVec_def = 0;     % Use default delta_pVec

notFit_def = [sigma_ellipseData, sigmacriterion_def, lapse_def, gamma_def, ...
    alpha_def, w_def, Tmax_def, pVec_def, beta_def, delta_def, delta_def, ...
    hRate_def, nu_p_def, delta_Tmax_def, delta_pVec_def];

inputParams(I_notFit) = notFit_def(I_notFit);

%inputParams

%% Compute the min negative log likelihood using BADS

% DEFINE PARAMETER BOUNDS
LB = paramBounds(:,1);
UB = paramBounds(:,2);
  
% Initialize parameters
params_initial = (LB+UB)/2;

% Run fit
[fitParams, ~] = fminsearchbnd(@(params2fit)changeprob_nll(params2fit, I_params, inputParams, NumTrials, mu, ...
    sigma_s, C, S, p_true, resp_obs, task, score, model, X), params_initial, LB, UB);

%% Recompute additional outputs from best fit params

inputParams(I_params) = fitParams;
% sigma = sqrt(sigma_s^2 + inputParams(1)^2);
[nLL,~,resp_model,p_estimate,~] = ...
    changeprob_nll(fitParams, I_params, inputParams, NumTrials, mu, sigma_s, C, S, p_true, resp_obs, task, score, model, X);


end

