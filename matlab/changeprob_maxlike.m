function [nLL, fitParams, resp_model, resp_obs, p_true, p_estimate, lmL, vbmc_fit] = changeprob_maxlike(model,...
    data, task, parameters, paramBounds, Nopts)
%% CHANGEPROB_MAXLIKE Max log likelihood for specified model

% [Documentation to be written]

if nargin < 2; data = []; end
if nargin < 3; task = []; end
if nargin < 4; parameters = []; end
if nargin < 5; paramBounds = []; end
if nargin < 6 || isempty(Nopts); Nopts = 20; end

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
    paramBounds_def = [1,30; 1,30; 0,0.1; -Inf,Inf; 0,1; 0,1; 2,200; ...
        1e-6,0.5-1e-6; 0,10; 1.01,5; 1.01,14; 0,1; 0,5; 1,200; 1e-6,0.5-1e-6];    
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
PLB = LB + 0.1*(UB-LB);
PUB = LB + 0.9*(UB-LB);
  
% Define function to fit
nLL_fun = @(params2fit)changeprob_nll(params2fit(:)', I_params, inputParams, NumTrials, mu, ...
    sigma_s, C, S, p_true, resp_obs, task, score, model, X);

badopts = bads('defaults');             % Get default parameters for BADS
badopts.UncertaintyHandling = 'no';     % Deterministic objective function

fitParams_list = NaN(Nopts,numel(PLB));
nLL_list = NaN(1,Nopts);

% Loop over different starting points
for iOpt = 1:Nopts
    % Random starting point
    params_initial = PLB + rand(size(PLB)).*(PUB - PLB);

    % Run fit
    [fitParams_list(iOpt,:),nLL_list(iOpt)] = bads(nLL_fun, params_initial(:)', LB(:)', UB(:)', PLB(:)', PUB(:)', [], badopts);
end

fitParams_list
nLL_list

[~,idx_best] = min(nLL_list);
fitParams = fitParams_list(idx_best,:);
    
%% Recompute additional outputs from best fit params

inputParams(I_params) = fitParams;
% sigma = sqrt(sigma_s^2 + inputParams(1)^2);
[nLL,~,resp_model,p_estimate,~] = nLL_fun(fitParams);

%% Compute estimate of marginal likelihood with VBMC

% Set VBMC options (these will probably become defaults in the near future)
vbmc_opts = vbmc('defaults');
% vbmc_opts.Plot = 'on';

vbmc_opts.NSgpMaxMain = 0;
vbmc_opts.SGDStepSize = 0.005;
vbmc_opts.RetryMaxFunEvals = vbmc_opts.MaxFunEvals;     % Retry variational optimization if first fails
vbmc_opts.gpMeanFun = 'negquadse';

% Choose starting point equal to maximum-likelihood fit, but ensure it is
% well inside bounds
LB_eff = LB + (PLB-LB)*0.01;
UB_eff = UB - (UB-PUB)*0.01;
startPoint = min(max(fitParams(:),LB_eff),UB_eff);

Nopts = min(max(ceil(Nopts/2),2),10);

% logprior = @(x) -sum(log(UB - LB));  % Uniform prior

% Uniform prior with slight smoothing very close to the edges (needed by VBMC)
logprior = @(x) log(msplinetrapezpdf(x,LB(:)',LB_eff(:)',UB_eff(:)',UB(:)'));

for iOpt = 1:Nopts
    [vp,~,~,~,output] = ...
        vbmc(@(params2fit) -nLL_fun(params2fit(:)')+logprior(params2fit(:)'), startPoint(:)', ...
            LB(:)', UB(:)', PLB(:)', PUB(:)', vbmc_opts);
        
    % Store results
    vbmc_fit.vps{iOpt} = vp;
    vbmc_fit.outputs{iOpt} = output;
end

[exitflag,best,idx_best,stats] = vbmc_diagnostics(vbmc_fit.vps,vbmc_fit.outputs);
vbmc_fit.diagnostics.exitflag = exitflag;
vbmc_fit.diagnostics.stats = stats;
vbmc_fit.diagnostics.best = best;
vbmc_fit.diagnostics.idx_best = idx_best;

if ~isempty(best)
    lmL = best.elbo; % Save estimate of log marginal likelihood
else
    lmL = NaN;       % Run failed!
end

end

