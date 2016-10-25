function [model_resp,log_P] = RLobserver_covert(parameters,sigma,dmu,X,p_initial,resp_obs,score)
%RLOBSERVER_COVERT Responses and log likelihoods for covert RL observer.
%
% ================ INPUT VARIABLES ====================
% PARAMETERS: sensory noise and learning rate. [1,2] (double)
% SIGMA: combined category and sensory noise. [scalar] (double)
% DMU: distance between category means. [scalar] (double)
% X: matrix of noisy measurements. [Nt,Ns] (double)
% P_INITIAL: initial probability. [scalar] (double)
% RESP_OBS: subject's categorization responses. [Nt,1] (double)
% SCORE: Trial feedback. [Nt,1] (double)
% 
% ================ OUTPUT VARIABLES ==================
% MODEL_RESP: model categorization responses, per trial/sample. [Nt,Ns] (double)
% LOG_PZ: log likelihood, per trial/sample. [Nt,Ns] (double)

[Nt,Ns] = size(X); % # trials and # samples

z_model = zeros(Nt,Ns);
model_resp = z_model;
log_P = zeros(Nt,Ns);

lambda = 1e-4; % Minimum lapse to avoid numerical trouble

for qq = 1:Ns
    x = X(:,qq); % vector of noisy measurements
    Gamma(1) = p_initial/(1-p_initial);
    z_model(1,qq) = sigma^2 * log(Gamma(1)) / dmu;
    if x(1) <= z_model(1,qq)
        model_resp(1,qq) = 1;
    else
        model_resp(1,qq) = 0;
    end
    model_resp(1,qq) = lambda/2 + (1-lambda)*model_resp(1,qq);
    log_P(1,qq) = log(model_resp(1,qq)).*(Chat(1)==1) + log(1-model_resp(1,qq)).*(Chat(1)~=1);
    for t = 2:Nt
        if score(t-1) == 0
            z_model(t,qq) = z_model(t-1,qq) + parameters(2)*(x(t)-z_model(t-1,qq));
        else
            z_model(t,qq) = z_model(t-1,qq);
        end
        if x(t) <= z_model(t,qq)
            model_resp(t,qq) = 1;
        else
            model_resp(t,qq) = 0;
        end
        model_resp(t,qq) = lambda/2 + (1-lambda)*model_resp(t,qq);
        log_P(t,qq) = log(model_resp(t,qq)).*(Chat(t)==1) + log(1-model_resp(t,qq)).*(Chat(t)~=1);
    end
end