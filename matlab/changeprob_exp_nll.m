function [nLL, rmse, p_estimate, resp_model] = changeprob_exp_nll(parameters, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task)
%CHANGEPROB_EXP_NLL Exponentially weight-moving average observer.
% (Documentation to be written.)
%
% Author:   Elyse Norton
% Email:    elyse.norton@gmail.com
% Date:     Oct/20/2016

% Parameter vector:
% #1 is SIGMA_ELLIPSE, #2 is SIGMA_CRITERION, #3 is LAPSE, #4 is GAMMA,
% #5 is ALPHA, and #6 is W
if nargin < 1; parameters = []; ...
        [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs] = changeprob_getSessionParameters(); task = 1; end

% Data struct or random seed for fake data generation
if nargin < 2; error('You must specify the session parameters.'); end

% Task (1 overt, 2 covert)
if nargin < 9 || isempty(task); task = 1; end
if task ~= 1 && task ~= 2; error('TASK can only be 1 (overt-criterion) or 2 (covert-criterion).'); end

%% Define parameters

switch task
    case 1  % Overt-criterion task    
        z_resp = resp_obs;   % Reported criterion
    case 2  % Covert-criterion task
        Chat = resp_obs;
end

% Observer model parameters

switch numel(parameters)
    case 0
        sigma_criterion = 5; % Default sigma_criterion
        lambda = 0; % Default lapse rate
        alpha = .2; % Defalt smoothing factor
        w = 1; % Default weight on probability estimate
    case 6
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = parameters(3);
        alpha = parameters(5);
        w = parameters(6);
    otherwise
        error('Parameters is a vector with exactly 0 or 6 inputs.');
end

%% Start loop over trials

p_initial = .5;
p_bias = .5; % Conservative bias towards a probability of .5 (i.e., neutral criterion)
z_model = zeros(1,NumTrials);
PChatA = z_model;
log_P = z_model;
beta = z_model;
p_estimate = zeros(NumTrials,1);
p_conservative = p_estimate;
p_estimate(1,:) = p_initial;
p_conservative(1,:) = w*p_estimate(1,:) + (1-w)*p_bias;
beta(1) = p_conservative(1,:)/(1-p_conservative(1,:));
z_model(:,1) = sigma^2 * log(beta(1)) / diff(mu);
if task == 1
    log_P(:,1) = -0.5*log(2*pi*sigma_criterion) - 0.5*((z_resp(1)-z_model(:,1))./sigma_criterion).^2;
else
    % PChatA(:,1) = 1 - normcdf(S(1), z_model(:,1), sigma_ellipse);
    PChatA(:,1) = 1 - 0.5*erfc( -(S(1) - z_model(1)) / (sqrt(2)*sigma_ellipse) );   % Faster implementation
    log_P(:,1) = log(PChatA(:,1)).*(Chat(1)==1) + log(1-PChatA(:,1)).*(Chat(1)==2);
end
Cprev = C(1);

for t = 2:NumTrials
    p_estimate(t,:) = p_estimate(t-1,:) + alpha*(Cprev-p_estimate(t-1,:));
    p_conservative(t,:) = w*p_estimate(t,:) + (1-w)*p_bias;
    beta(t) = p_conservative(t,:)/(1-p_conservative(t,:));
    z_model(t) = sigma^2 * log(beta(t)) / diff(mu);
    if task == 1
        log_P(:,t) = -0.5*log(2*pi*sigma_criterion) - 0.5*((z_resp(t)-z_model(:,t))./sigma_criterion).^2;
        if lambda > 0
            log_P(:,t) = log(lambda/360 + (1-lambda)*exp(log_P(:,t)));
        end
    else
        % PChatA(:,t) = 1-normcdf(S(t), z_model(t), sigma_ellipse);
        PChatA(:,t) = 1 - 0.5*erfc( -(S(t) - z_model(t)) / (sqrt(2)*sigma_ellipse) );   % Faster implementation
        if lambda > 0
            PChatA(:,t) = lambda/2 + (1-lambda)*PChatA(:,t);
        end
        log_P(:,t) = log(PChatA(:,t)).*(Chat(t)==1) + log(1-PChatA(:,t)).*(Chat(t)~=1);
    end
    Cprev = C(t);
end

if task == 1
    resp_model = z_model;
else
    resp_model = PChatA;
end

nLL = -nansum(log_P);
rmse = sqrt(mean((p_estimate - p_true).^2));

end

