function [nLL, rmse, p_estimate, p_true, resp_model, resp_obs] = changeprob_RLProb_nll(parameters, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task)
%CHANGEPROB_RLProb_NLL Reinforcement-learning model on probability.
% (Documentation to be written.)
%
% Author:   Elyse Norton
% Email:    elyse.norton@gmail.com
% Date:     Dec/12/2016

% Generate fake data and set default parameter vector
% Parameter vector: #1 is SIGMA_ELLIPSE, #2 is SIGMA_CRITERION, #3 is LAPSE, 
% #4 is GAMMA, #5 is ALPHA, #6 is W, #7 is Tmax
if nargin < 1 || isempty(parameters)
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs] = changeprob_getSessionParameters();
    parameters = [sigma_ellipse, sigma_ellipse, 1e-4, Inf, .2, 1, 0];
    task = 1; 
end

% Generate fake data using parameters specified
if nargin < 2
    task = 1;
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs] = changeprob_getSessionParameters([], task, [parameters(1), parameters(2)]);
end

% Make sure all the experiment variables are specified
if nargin > 1 && nargin < 9
    error('Missing input arguments. You must specify all or no experimental variables.');
end

% Task (1 overt, 2 covert)
if nargin == 9 || isempty(task); task = 1; end
if task ~= 1 && task ~= 2 && task ~=3; error('TASK can only be 1 (overt-criterion), 2 (covert-criterion), or 3 (mixed design).'); end

% Observer model parameters

sigma_ellipse = parameters(1);
sigma_criterion = parameters(2);
lambda = parameters(3);
alpha = parameters(5);
w = parameters(6);

%% Start loop over trials

p_initial = .5;
p_bias = .5; % Conservative bias towards a probability of .5
z_model = zeros(1,NumTrials);
PChatA = ones(1,NumTrials);
log_P = zeros(1,NumTrials);
beta = zeros(1,NumTrials);
p_estimate = zeros(NumTrials,1);
p_conservative = zeros(NumTrials,1);
p_estimate(1,:) = p_initial;
p_conservative(1,:) = w*p_estimate(1,:) + (1-w)*p_bias;
beta(1) = p_conservative(1,:)/(1-p_conservative(1,:));
z_model(:,1) = sigma^2 * log(beta(1)) / diff(mu);
if task == 1
    log_P(:,1) = -0.5*log(2*pi*sigma_criterion) - 0.5*((resp_obs(1)-z_model(:,1))./sigma_criterion).^2;
else
    % PChatA(:,1) = 1 - normcdf(S(1), z_model(:,1), sigma_ellipse);
    PChatA(:,1) = 1 - 0.5*erfc( -(S(1) - z_model(1)) / (sqrt(2)*sigma_ellipse) );   % Faster implementation            
    log_P(:,1) = log(PChatA(:,1)).*(resp_obs(1)==1) + log(1-PChatA(:,1)).*(resp_obs(1)~=1);
end
Cprev = C(1);

for t = 2:NumTrials
    if score(t-1) == 0
        p_estimate(t,:) = p_estimate(t-1,:) + alpha*(Cprev-p_estimate(t-1,:));
    else
        p_estimate(t,:) = p_estimate(t-1,:);
    end
    p_conservative(t,:) = w*p_estimate(t,:) + (1-w)*p_bias;
    beta(t) = p_conservative(t,:)/(1-p_conservative(t,:));
    z_model(t) = sigma^2 * log(beta(t)) / diff(mu);
    if task == 1 || and(task == 3, mod(t,5) == 0)
        log_P(:,t) = -0.5*log(2*pi*sigma_criterion) - 0.5*((resp_obs(t)-z_model(:,t))./sigma_criterion).^2;
        if lambda > 0
            log_P(:,t) = log(lambda/360 + (1-lambda)*exp(log_P(:,t)));
        end
    else
        % PChatA(:,t) = 1-normcdf(S(t), z_model(t), sigma_ellipse);
        PChatA(:,t) = 1 - 0.5*erfc( -(S(t) - z_model(t)) / (sqrt(2)*sigma_ellipse) );   % Faster implementation            
        if lambda > 0
            PChatA(:,t) = lambda/2 + (1-lambda)*PChatA(:,t);
        end
        log_P(:,t) = log(PChatA(:,t)).*(resp_obs(t)==1) + log(1-PChatA(:,t)).*(resp_obs(t)~=1);
    end
    Cprev = C(t);
end

if task == 1
    resp_model = z_model;
elseif task ==2
    resp_model = PChatA;
else
    resp_model = PChatA;
    resp_model(5:5:NumTrials) = z_model(5:5:NumTrials);
end

nLL = -nansum(log_P);
rmse = sqrt(mean((p_estimate - p_true).^2));

end