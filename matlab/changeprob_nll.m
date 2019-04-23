function [nLL, rmse, resp_model, p_estimate, post, v_estimate] = changeprob_nll(params, idx_params, ...
    inputParams, NumTrials, mu, sigma_s, C, S, p_true, resp_obs, task, score, model, X)
%CHNAGEPROB_NLL Computes the negative log likelihood for the specified
%model
%(Documentation to be written)

inputParams(idx_params) = params;
sigma = sqrt(sigma_s^2 + inputParams(1)^2);
v_estimate = [];

switch model
    case 1
        if inputParams(7) == 0
            prior_rl = [];
        elseif (inputParams(7) ~= 0) && (inputParams(14) == 0)
            prior_rl = [max(1,inputParams(7)*2/3),inputParams(7)];
        else
            prior_rl = [(inputParams(7)-1), inputParams(7)-1+inputParams(14)];
        end
        if inputParams(8) == 0
            p_vec = [];
        elseif (inputParams(8) == 1) && (inputParams(15) == 0)
            p_vec = linspace(inputParams(8), 1-inputParams(8), 5);
        else
            p_vec = linspace(inputParams(8), 0.5+inputParams(15), 5);
        end
        if inputParams(9) == 0
            beta_hyp = [];
        else
            beta_hyp = inputParams(9)^2;
        end
        [nLL, rmse, p_estimate, resp_model, post] = ChangeProb_bocpd_nll_v2(inputParams(1:6), NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task, prior_rl, p_vec, beta_hyp);
    case 2
        [nLL, rmse, p_estimate, resp_model] = changeprob_fixed_nll(inputParams(1:6), NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task);
        post = [];
    case 3
        [nLL, rmse, p_estimate, resp_model] = changeprob_exp_nll(inputParams(1:6), NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task);
        post = [];
    case 4
        [nLL, rmse, p_estimate, resp_model] = changeprob_RLProb_nll(inputParams(1:6), NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task);
        post = [];
    case 5
        switch task
            case 1
                if isempty(X)
                    X = bsxfun(@plus, S, inputParams(1)*randn(numel(S), 5000));
                end
                [resp_model, logP] = RLobserver_overt_mex([inputParams(5), inputParams(2)],sigma,diff(mu),X,.5,resp_obs,score);
            case 2
                X = bsxfun(@plus, S, inputParams(1)*randn(numel(S), 5000));
                [resp_model, logP] = RLobserver_covert_mex([inputParams(1), inputParams(5)],sigma,diff(mu),X,.5,resp_obs,score);
            case 3
                [resp_model, logP] = RLobserver_mixed([inputParams(2), inputParams(5)],sigma,diff(mu),X,.5,resp_obs,score);
        end
        % Marginalize over x (in likelihood space NOT log likelihood space)
        maxLL = max(logP(:));
        likelihood = mean(exp(logP-maxLL),2); % p(resp | x, theta)
        nLL = -nansum(log(likelihood)+maxLL);
        resp_model = mean(resp_model,2)';
        rmse = [];
        p_estimate = [];
        post = [];
    case 6
        [nLL, rmse, p_estimate, resp_model] = changeprob_gold_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task);
        post = [];
    case 7
        [nLL, rmse, p_estimate, resp_model, v_estimate] = changeprob_behrens(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task);
        post = [];
    case 8
        [nLL, rmse, p_estimate, resp_model, v_estimate] = changeprob_behrens_jump(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task);
        post = [];
    otherwise
        error('Cannot compute the negative log likelihood for the chosen model.');
end

end
