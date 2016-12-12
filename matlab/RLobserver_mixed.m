function [model_resp,log_P] = RLobserver_mixed(parameters,sigma,dmu,X,p_initial,resp_obs,score)
%RLOBSERVER_COVERT Responses and log likelihoods for covert RL observer.
%
% ================ INPUT VARIABLES ====================
% PARAMETERS: adjustment noise and learning rate. [1,2] (double)
% SIGMA: combined category and sensory noise. [scalar] (double)
% DMU: distance between category means. [scalar] (double)
% X: matrix of noisy measurements. [Nt,Ns] (double)
% P_INITIAL: initial probability. [scalar] (double)
% RESP_OBS: subject's responses. [Nt,1] (double)
% SCORE: Trial feedback. [Nt,1] (double)
% 
% ================ OUTPUT VARIABLES ==================
% MODEL_RESP: model responses, per trial/sample. [Nt,Ns] (double)
% LOG_P: log likelihood, per trial/sample. [Nt,Ns] (double)

[Nt,Ns] = size(X); % # trials and # samples

z_model = zeros(Nt,Ns);
model_resp = z_model;
log_P = zeros(Nt,Ns);

lambda = 1e-2; % Minimum lapse to avoid numerical trouble (CHANGED)

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
    log_P(1,qq) = log(model_resp(1,qq)).*(resp_obs(1)==1) + log(1-model_resp(1,qq)).*(resp_obs(1)~=1);
    xprev = x(1);
    for t = 2:Nt
        if score(t-1) == 0
            % This version is 25 times slower
            % z_model(t,qq) = z_model(t-1,qq) + parameters(2)*(x(t-1)-z_model(t-1,qq)); 
            z_model(t,qq) = z_model(t-1,qq) + parameters(2)*(xprev-z_model(t-1,qq));
        else
            z_model(t,qq) = z_model(t-1,qq);
        end
        if mod(t,5) == 0
            model_resp(t,qq) = z_model(t,qq);
            log_P(t,qq) = -0.5*log(2*pi*parameters(1)) - 0.5*((resp_obs(t)-z_model(t,qq))./parameters(1)).^2;
            if lambda > 0
                log_P(t,qq) = log(lambda/360 + (1-lambda)*exp(log_P(t,qq)));
            end
        else  
            if x(t) <= z_model(t,qq)
                model_resp(t,qq) = 1;
            else
                model_resp(t,qq) = 0;
            end
            model_resp(t,qq) = lambda/2 + (1-lambda)*model_resp(t,qq);
            log_P(t,qq) = log(model_resp(t,qq)).*(resp_obs(t)==1) + log(1-model_resp(t,qq)).*(resp_obs(t)~=1);
        end
        xprev = x(t);
    end
end
