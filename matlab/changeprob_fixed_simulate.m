function dataSim = changeprob_fixed_simulate(data, task, model, parameters)
%CHANGEPROB_FIXED_SIMULATE Simulates a fixed criterion observer.
% (Documentation to be written.)
%
% Author:   Elyse Norton
% Email:    elyse.norton@gmail.com
% Date:     Dec/30/2016

% Check input arguments

% Check for data file
if nargin < 1; data = []; task = 1; parameters = []; end

% Determine task to simulate
if nargin < 2; task = 1; parameters = []; end
if isempty(task); task = 1; end
if task ~= 1 && task ~= 2 && task ~=3; error('TASK can only be 1 (overt-criterion), 2 (covert-criterion), or 3 (mixed design).'); end

if nargin < 3; model = {'fixed'}; parameters = []; end

% Parameter vector: SIGMA_ELLIPSE or SIGMA_CRITERION
if nargin < 4 || isempty(parameters); parameters = 5; end

% Get session parameters
if isempty(data) 
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs] = changeprob_getSessionParameters([], task, [parameters(1), parameters(1)]); % Generate fake data using specified parameters
else
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs] = changeprob_getSessionParameters(data, task); % Determine session parameters from provided data set
end

switch numel(parameters)
    case 1
        if task ~= 2
            sigma_criterion = parameters(1);
        else
            sigma_ellipse = parameters(1);
        end
    otherwise
        error('Too many input parameters');
end      

%% Start loop over trials

p_initial = .5;
p_estimate = bsxfun(@plus, zeros(NumTrials,1), p_initial); % Assume probability is fixed throughout the experimental block
beta = p_estimate./(1-p_estimate);
z_model = sigma^2 * log(beta)/diff(mu);
resp_obs = zeros(NumTrials,1);
score = resp_obs;

switch task
    case 1
        resp_obs = z_model + randn(NumTrials,1)*sigma_criterion;
        score = double(or(C == 1 & S <= resp_obs, C == 0 & S > resp_obs));
    case 2
        X = S + randn(NumTrials,1)*sigma_ellipse;
        resp_obs = double(X <= z_model);
        score = double(C == resp_obs);
    case 3
        idx_Overt = 5:5:NumTrials;
        resp_obs(idx_Overt) = z_model(idx_Overt) + randn(numel(idx_Overt),1)*sigma_criterion;
        score(idx_Overt) = double(or(C(idx_Overt) == 1 & S(idx_Overt) <= resp_obs(idx_Overt), C(idx_Overt) == 0 & S(idx_Overt) > resp_obs(idx_Overt)));
        idx_Covert = 1:NumTrials;
        idx_Covert(idx_Overt) = [];
        X = S(idx_Covert) + randn(numel(idx_Covert),1)*sigma_ellipse;
        resp_obs(idx_Covert) = double(X <= z_model(idx_Covert));
        score(idx_Covert) = double(C(idx_Covert) == resp_obs(idx_Covert));
end

dataSim.NumTrials = NumTrials;
dataSim.mu = mu;
dataSim.sigma = sigma;
dataSim.sigmaEllipse = sigma_ellipse;
dataSim.Category = C;
dataSim.Stimulus = S;
dataSim.pA = p_true;
dataSim.response = resp_obs;
dataSim.score = score;

end