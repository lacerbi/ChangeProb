function [nLL, rmse, p_estimate, resp_model] = changeprob_fixed_nll(parameters, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task)
%CHANGEPROB_FIXED_NLL Fixed criterion model.
% (Documentation to be written.)
%
% Author:   Elyse Norton
% Email:    elyse.norton@gmail.com
% Date:     Dec/12/2016

% Generate fake data and set default parameter vector
% Parameter vector: #1 is SIGMA_ELLIPSE, #2 is SIGMA_CRITERION, #3 is LAPSE, 
% #4 is GAMMA, #5 is ALPHA, and #6 is W
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

%% Start loop over trials

p_initial = .5;
p_estimate = bsxfun(@plus, zeros(NumTrials,1), p_initial); % Assume probability is fixed throughout the experimental block
beta = p_estimate./(1-p_estimate);
z_model = sigma^2 * log(beta)/diff(mu);
if task == 3
    idx_Overt = 5:5:NumTrials;
    idx_Covert = 1:NumTrials;
    idx_Covert(idx_Overt) = [];
end

switch task
    case 1
        log_P = -0.5*log(2*pi*sigma_criterion) - 0.5*((resp_obs-z_model)./sigma_criterion).^2;
        if lambda > 0
            log_P = log(lambda/360 + (1-lambda).*exp(log_P));
        end
        resp_model = z_model;
    case 2
        PChatA = 1 - normcdf(S, z_model, sigma_ellipse);
        if lambda > 0
            PChatA = lambda/2 + (1-lambda).*PChatA;
        end
        log_P = log(PChatA).*(resp_obs==1) + log(1-PChatA).*(resp_obs~=1);
        resp_model = PChatA;
    case 3
        log_P = zeros(1, NumTrials);
        resp_model = zeros(1, NumTrials);
        
        % Log probability of overt responses
        log_P(idx_Overt) = -0.5*log(2*pi*sigma_criterion) - 0.5*((resp_obs(idx_Overt)-z_model(idx_Overt))./sigma_criterion).^2;
        if lambda > 0
            log_P(idx_Overt) = log(lambda/360 + (1-lambda).*exp(log_P(idx_Overt)));
        end
        
        % Log probability of covert responses
        PChatA = 1 - normcdf(S(idx_Covert), z_model(idx_Covert), sigma_ellipse);
        if lambda > 0
            PChatA = lambda/2 + (1-lambda).*PChatA;
        end
        log_P(idx_Covert) = log(PChatA).*(resp_obs(idx_Covert)==1) + log(1-PChatA).*(resp_obs(idx_Covert)~=1);
        
        resp_model(idx_Covert) = PChatA;
        resp_model(idx_Overt) = z_model(idx_Overt);
end

nLL = -nansum(log_P);
rmse = sqrt(mean((p_estimate - p_true).^2));

end