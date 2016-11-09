function [nLL, rmse, resp_model, p_estimate, post] = changeprob_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, task, score, model, X)
%CHNAGEPROB_NLL Computes the negative log likelihood for the specified
%model
%(Documentation to be written)

switch model
    case 1
        if inputParams(7) == 0
            prior_rl = [];
        else
            prior_rl = [max(1,floor(inputParams(7)*2/3)),inputParams(7)];
        end
        [nLL, rmse, p_estimate, resp_model, post] = ChangeProb_bocpd_nll_v2(inputParams(1:4), NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task, prior_rl);
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
                [resp_model, logP] = RLobserver_overt_mex([inputParams(5), inputParams(2)],sigma,diff(mu),X,.5,resp_obs,score);
            case 2
                X = bsxfun(@plus, S, inputParams(1)*randn(numel(S), 5000));
                [resp_model, logP] = RLobserver_covert_mex([inputParams(1), inputParams(5)],sigma,diff(mu),X,.5,resp_obs,score);
        end
        % Marginalize over x (in likelihood space NOT log likelihood space)
        maxLL = max(logP(:));
        likelihood = mean(exp(logP-maxLL),2); % p(resp | x, theta)
        nLL = -nansum(log(likelihood)+maxLL);
        resp_model = mean(resp_model,2)';
        rmse = [];
        p_estimate = [];
        post = [];
    otherwise
        error('Cannot compute the negative log likelihood for the chosen model.');
end

end
