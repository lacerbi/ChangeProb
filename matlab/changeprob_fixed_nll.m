function [nLL, rmse, p_estimate, resp_model] = changeprob_fixed_nll(parameters, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task)
%CHANGEPROB_FIXED_NLL Fixed criterion model.
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
    case 6
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = parameters(3);
    otherwise
        error('Parameters is a vector with exactly 0 or 6 inputs.');
end

%% Start loop over trials

p_initial = .5;
p_estimate = zeros(NumTrials,1) + p_initial; % Assume probability is fixed throughout the experimental block
beta = p_estimate./(1-p_estimate);
z_model = sigma^2 * log(beta)/diff(mu);

switch task
    case 1
        log_P = -0.5*log(2*pi*sigma_criterion) - 0.5*((z_resp-z_model)./sigma_criterion).^2;
        if lambda > 0
            log_P = log(lambda/360 + (1-lambda).*exp(log_P));
        end
        resp_model = z_model;
    case 2
        PChatA = 1 - normcdf(S, z_model, sigma_ellipse);
        log_P = log(PChatA).*(Chat==1) + log(1-PChatA).*(Chat==2);
        if lambda > 0
                PChatA = lambda/2 + (1-lambda).*PChatA;
        end
        resp_model = PChatA;
end

nLL = -nansum(log_P);
rmse = sqrt(mean((p_estimate - p_true).^2));

end