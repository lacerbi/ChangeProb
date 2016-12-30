function dataSim = changeprob_RLcriterion_simulate(data, task, parameters)
%CHANGEPROB_EXP_SIMULATE Simulates an RL criterion observer.
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

% Parameter vector: #1 is SIGMA_ELLIPSE / SIGMA_CRITERION, #2 is ALPHA
if nargin < 3 || isempty(parameters); parameters = [5 .4]; end

% Get session parameters
if isempty(data) 
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs] = changeprob_getSessionParameters([], task, [parameters(1), parameters(1)]); % Generate fake data using specified parameters
else
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs] = changeprob_getSessionParameters(data, task); % Determine session parameters from provided data set
end

switch numel(parameters)
    case 1
        error('Not enough input parameters');
    case 2
        if task ~= 2
            sigma_criterion = parameters(1);
        else
            sigma_ellipse = parameters(1);
        end
        alpha = parameters(2);
    otherwise
        error('Too many input parameters');
end     

%% Start loop over trials

p_initial = .5;
z_model = zeros(1,NumTrials);
resp_obs = z_model;
score = z_model;
beta = p_initial/(1-p_initial);
z_model(:,1) = sigma^2 * log(beta) / diff(mu);
X = S + randn(NumTrials,1)*sigma_ellipse;

if task == 1
    resp_obs(:,1) = z_model(:,1) + randn(1,1)*sigma_criterion;
    if and(C(1) == 1, S(1) <= resp_obs(:,1)) || and(C(1) == 0, S(1) > resp_obs(:,1))
        score(:,1) = 1;
    else
        score(:,1) = 0;
    end
else
    if X(1) <= z_model(:,1)
        resp_obs(:,1) = 1;
    else
        resp_obs(:,1) = 0;
    end
    score(:,1) = double(resp_obs(:,1) == C(1));
end

Xprev = X(1);

for t = 2:NumTrials
    if score(:,t-1) == 0
        z_model(:,t) = z_model(:,t-1) + alpha*(Xprev-z_model(:,t-1));
    else
        z_model(:,t) = z_model(:,t-1);
    end
    if task == 1 || and(task == 3, mod(t,5)==0)
        resp_obs(:,t) = z_model(:,t) + randn(1,1)*sigma_criterion;
        if and(C(t) == 1, S(t) <= resp_obs(:,t)) || and(C(t) == 0, S(t) > resp_obs(:,t))
            score(:,t) = 1;
        else
            score(:,t) = 0;
        end
    else
        if X(t) <= z_model(:,t)
            resp_obs(:,t) = 1;
        end
        score(:,t) = double(resp_obs(:,t) == C(t));
    end
    Xprev = X(t);
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