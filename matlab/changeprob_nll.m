function [nLL, rmse, resp_model, p_estimate, post] = changeprob_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, task, score, model)
%CHNAGEPROB_NLL Computes the negative log likelihood for the specified
%model
%(Documentation to be written)

switch model
    case 1
        [nLL, rmse, p_estimate, resp_model, post] = ChangeProb_bocpd_nll_v2(inputParams(1:4), NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task);
    case 2
        [nLL, rmse, p_estimate, resp_model] = changeprob_fixed_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task);
        post = [];
    case 3
        [nLL, rmse, p_estimate, resp_model] = changeprob_exp_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task);
        post = [];
    case 4
        [nLL, rmse, p_estimate, resp_model] = changeprob_RLProb_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task);
        post = [];
    case 5
        switch task
            case 1
                X = S + inputParams(1)*randn(numel(S), 1000);
                [resp_model, logP] = RLobserver_overt_mex(inputParams,sigma,diff(mu),X,.5,resp_obs,score);
            case 2
                X = S + inputParams(1)*randn(numel(S), 1000);
                [resp_model, logP] = RLobserver_covert_mex(inputParams,sigma,diff(mu),X,.5,resp_obs,score);
        end
        % Marginalize over x (in likelihood space NOT log likelihood space)
        logLikelihood = nansum(logP); % Log likelihood for each x vector
        likelihood = mean(exp((logLikelihood-max(logLikelihood)))); % p(resp | x, theta)
        rmse = [];
        p_estimate = [];
        post = [];
    otherwise
        error('Cannot compute the negative log likelihood for the chosen model.');
end

end
