function [model_resp,log_Pz] = RLobserver_overt(parameters,sigma,dmu,X,p_initial,z_resp,score)
%RLOBSERVER_OVERT Responses and log likelihoods for overt RL observer.
%
% ================ INPUT VARIABLES ====================
% PARAMETERS: learning rate and adjustment noise. [1,2] (double)
% SIGMA: combined category and sensory noise. [scalar] (double)
% DMU: distance between category means. [scalar] (double)
% X: matrix of noisy measurements. [Nt,Ns] (double)
% P_INITIAL: initial probability. [scalar] (double)
% Z_RESP: subject's criterion responses. [Nt,1] (double)
% SCORE: Trial feedback. [Nt,1] (double)
% 
% ================ OUTPUT VARIABLES ==================
% MODEL_RESP: model criterion responses, per trial/sample. [Nt,Ns] (double)
% LOG_PZ: log likelihood, per trial/sample. [Nt,Ns] (double)

[Nt,Ns] = size(X); % # trials and # samples

model_resp = zeros(Nt,Ns);
log_Pz = zeros(Nt,Ns);

for qq = 1:Ns
    x = X(:,qq); % vector of noisy measurements
    Gamma(1) = p_initial/(1-p_initial);
    model_resp(1,qq) = sigma^2 * log(Gamma(1)) / dmu;
    log_Pz(1,qq) = -0.5*log(2*pi*parameters(2)) - 0.5*((z_resp(1)-model_resp(1,qq))./parameters(2)).^2;
    xprev = x(1);
    for t = 2:Nt
        if score(t-1) == 0
            model_resp(t,qq) = model_resp(t-1,qq) + parameters(1)*(xprev-model_resp(t-1,qq));
        else
            model_resp(t,qq) = model_resp(t-1,qq);
        end
        xprev = x(t);
        log_Pz(t,qq) = -0.5*log(2*pi*parameters(2)) - 0.5*((z_resp(t)-model_resp(t,qq))./parameters(2)).^2;
    end
end
