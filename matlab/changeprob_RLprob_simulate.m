function dataSim = changeprob_RLprob_simulate(data, task, model, parameters)
%CHANGEPROB_EXP_SIMULATE Simulates an RL probability observer.
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

if nargin < 3; model = {'RL_probability'}; parameters = []; end

% Parameter vector: #1 is SIGMA_ELLIPSE / SIGMA_CRITERION, #2 is ALPHA, #3
% is W (optional)
if nargin < 4 || isempty(parameters); parameters = [5 .1 1]; end

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
        w = 1;
    case 3
        if task ~= 2
            sigma_criterion = parameters(1);
        else
            sigma_ellipse = parameters(1);
        end
        alpha = parameters(2);
        w = parameters(3);
    otherwise
        error('Too many input parameters');
end     

%% Start loop over trials

p_initial = .5;
p_bias = .5; % Conservative bias towards a probability of .5 (i.e., neutral criterion) when w < 1
z_model = zeros(1,NumTrials);
resp_obs = z_model;
score = z_model;
beta = z_model;
p_estimate = zeros(NumTrials,1);
p_conservative = p_estimate;
p_estimate(1,:) = p_initial;
p_conservative(1,:) = w*p_estimate(1,:) + (1-w)*p_bias;
beta(1) = p_conservative(1,:)/(1-p_conservative(1,:));
z_model(:,1) = sigma^2 * log(beta(1)) / diff(mu);

if task == 1
    resp_obs(:,1) = z_model(:,1) + randn(1,1)*sigma_criterion;
    if and(C(1) == 1, S(1) <= resp_obs(:,1)) || and(C(1) == 0, S(1) > resp_obs(:,1))
        score(:,1) = 1;
    else
        score(:,1) = 0;
    end
else
    if S(1)+randn(1,1)*sigma_ellipse <= z_model(:,1)
        resp_obs(:,1) = 1;
    else
        resp_obs(:,1) = 0;
    end
    score(:,1) = double(resp_obs(:,1) == C(1));
end

Cprev = C(1);

for t = 2:NumTrials
    if score(:,t-1) == 0
        p_estimate(t,:) = p_estimate(t-1,:) + alpha*(Cprev-p_estimate(t-1,:));
    else
        p_estimate(t,:) = p_estimate(t-1,:);
    end
    p_conservative(t,:) = w*p_estimate(t,:) + (1-w)*p_bias;
    beta(t) = p_conservative(t,:)/(1-p_conservative(t,:));
    z_model(t) = sigma^2 * log(beta(t)) / diff(mu);
    if task == 1 || and(task == 3, mod(t,5)==0)
        resp_obs(:,t) = z_model(:,t) + randn(1,1)*sigma_criterion;
        if and(C(t) == 1, S(t) <= resp_obs(:,t)) || and(C(t) == 0, S(t) > resp_obs(:,t))
            score(:,t) = 1;
        else
            score(:,t) = 0;
        end
    else
        if S(t)+randn(1,1)*sigma_ellipse <= z_model(:,t)
            resp_obs(:,t) = 1;
        end
        score(:,t) = double(resp_obs(:,t) == C(t));
    end
    Cprev = C(t);
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