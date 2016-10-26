function [marginalLikelihood, modelPost, nLL, rmse, fitParams, resp_model,...
    resp_obs, p_true, p_estimate] = changeprob_mL(model, data, task, parameters, gridSize, paramBounds)

%% CHANGEPROB_ML Computes marginal likelihood for the specified model

    % Here we compute the marginal likelihood of the specfied model. 
    % To calculate the marginal likelihood, we multiply the 
    % likelihood of the data given the model by the prior of the model and 
    % integrate (i.e., sum) over the model parameters.
    
    % INPUT:
        % model: a character string indicating the model you want to fit, 
        % which is limited to the following:
            % 'idealBayesian'
            % 'fixed'
            % 'exponential'
            % 'RL_probability'
            % 'RL_criterion'
        % data: data struct from changing probability experiment
        % task: lets you choose which task to fit
            % 1 - overt-criterion task
            % 2 - covert-criterion task
        % Vector indicating the parameters to-be-fit (1 - fit, 0 - not fit)
            % [sigma_ellipse, sigma_criterion, lapse, gamma, alpha, w]
        % gridSize: size of parameter grid (e.g., [n x m x p])
        % paramBounds: lower and upper parameter bounds (numel(gridSize) x 2) 
        
    % OUTPUT:
        % marginalLikelihood: a measure of model fit
        % modelPost: proportional to the p(parameters | data, model)
        % nLL: negative log likelihood that correponds to the best fitting
        % parameters
        % rmse: root mean squared error between the model probability and
        % the true probability (Note: does not apply to the RL_criterion model)
        % fitParams: best fitting model parameters computed by taking the
        % MAP of modelPost
        % resp_model: predicted criterion
        % resp_obs: observer's criterion
        % p_true: true probability values
        % p_estimate: model's estimate of probability (Note: does not apply to
        % the RL_criterion model)
        
    % Code adapted from Luigi Acerbi's changeprob_bocpd_nll.m 
    
    % Author:   Elyse norton
    % Email:    elyse.norton@gmail.com
    % Date:     10/20/2016
    tic
    % Model to be fit
    if nargin < 1; error('Please indicate the model you want to fit.'); end
    potentialModels = {'idealBayesian', 'fixed', 'exponential', 'RL_probability', ...
        'RL_criterion'};
    model = find(strcmp(model, potentialModels)); % recode model to numeric value
    
    % Data struct or random seed for fake data generation
    if nargin < 2 || isempty(data); data = 0; end
    if isnumeric(data); rng(data); data = []; end

    % Task (1 overt, 2 covert)
    if nargin < 3 || isempty(task); task = 1; end % overt-criterion task is the default
    if task ~= 1 && task ~= 2; error('TASK can only be 1 (overt) or 2 (covert).'); end

    % Parameters to be fit
    if nargin < 4 || isempty(parameters)
        switch task
            case 1
                if model < 3
                    parameters = [0 1 0 0 0 0];
                else
                    parameters = [0 1 0 0 1 0];
                end
            case 2
                if model < 3
                    parameters = [1 0 0 0 0 0];
                else
                    parameters = [1 0 0 0 1 0];
                end
        end
        if nargin < 4
            gridSize = [];
            paramBounds = [];
        end
    elseif sum(parameters) == 0
        error('You must specify at least one parameter to fit.')
    end
    
    paramNames = {'sigma_ellipse', 'sigma_criterion', 'lambda', 'gamma', 'alpha', 'w'};
    NumParams = sum(parameters);
    I_params = find(parameters ~= 0);
    % Sampling rate of parameter grid
    if nargin < 5 || isempty(gridSize)
        gridSize = zeros(1, NumParams) + 100; % Default grid is 100 x 100
        if nargin < 5
            paramBounds = [];
        end
    elseif numel(gridSize) ~= NumParams
        error('Matrix dimensions do not agree. numel(gridSize) must equal sum(parameters).');
    end
    
    % Lower and upper parameter bounds
    if nargin < 6 || isempty(paramBounds)
        paramBounds = zeros(NumParams, 2);
        for u = 1:NumParams
            if or(I_params(u) == 1, I_params(u) == 2)
                paramBounds(u,:) = [1, 30];
            elseif I_params(u) == 3
                paramBounds(u,:) = [0, .1];
            elseif I_params(u) == 4
                paramBounds(u,:) = [-Inf, Inf];
            else
                paramBounds(u,:) = [0, 1];
            end
        end
    end
    
    if NumParams ~= size(paramBounds,1)
        error('Please specify parameter bounds for each parameter to be fit.');
    end
    
    do_plot = nargout == 0; % If no outputs, make a plot
    
    %% Get session parameters
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs, score] = changeprob_getSessionParameters(data, task);
    
    %% Observer model parameters
    params2fit = zeros(NumParams, gridSize(1));
    for uu = 1:NumParams
        if or(I_params(uu) == 1, I_params(uu) == 2)
            params2fit(uu,:) = log(linspace(paramBounds(uu,1), paramBounds(uu,2), gridSize(uu))); % sigma_noise
        elseif I_params(uu) == 3
            params2fit(uu,:) = linspace(paramBounds(uu,1), paramBounds(uu,2), gridSize(uu)); % lambda
        elseif I_params(uu) == 4
            params2fit(uu,:) = linspace(paramBounds(uu,1), paramBounds(uu,2), gridSize(uu)); % gamma
        elseif I_params(uu) == 5
            params2fit(uu,:) = linspace(paramBounds(uu,1), paramBounds(uu,2), gridSize(uu)); % alpha
        else
            params2fit(uu,:) = linspace(paramBounds(uu,1), paramBounds(uu,2), gridSize(uu)); % w
        end
    end
    
    % Create parameter grid (max 5 parameters can be fit at once)
    Grid = cell(gridSize);
    switch NumParams
        case 1
            Grid = params2fit(1,:);
        case 2
            firstParam = params2fit(1,:);
            secondParam = params2fit(2,:);
            for i = 1:gridSize(1)
                for j = 1:gridSize(2)
                    Grid{i,j} = [firstParam(i), secondParam(j)];
                end
            end
        case 3
            firstParam = params2fit(1,:);
            secondParam = params2fit(2,:);
            thirdParam = params2fit(3,:);
            for i = 1:gridSize(1)
                for j = 1:gridSize(2)
                    for k = 1:gridSize(3)
                        Grid{i,j,k} = [firstParam(i), secondParam(j), thirdParam(k)];
                    end
                end
            end
        case 4
            firstParam = params2fit(1,:);
            secondParam = params2fit(2,:);
            thirdParam = params2fit(3,:);
            fourthParam = params2fit(4,:);
            for i = 1:gridSize(1)
                for j = 1:gridSize(2)
                    for k = 1:gridSize(3)
                        for q = 1:gridSize(4)
                            Grid{i,j,k,q} = [firstParam(i), secondParam(j), thirdParam(k), fourthParam(q)];
                        end
                    end
                end
            end
        case 5
            firstParam = params2fit(1,:);
            secondParam = params2fit(2,:);
            thirdParam = params2fit(3,:);
            fourthParam = params2fit(4,:);
            fifthParam = params2fit(5,:);
            for i = 1:gridSize(1)
                for j = 1:gridSize(2)
                    for k = 1:gridSize(3)
                        for q = 1:gridSize(4)
                            for p = gridSize(5)
                                Grid{i,j,k,q,p} = [firstParam(i), secondParam(j), thirdParam(k), fourthParam(q), fifthParam(p)];
                            end
                        end
                    end
                end
            end
        otherwise
            error('changeprob_mL can only fit up to 5 parameters');
    end
    
    % Default parameter values for those not fitted
    inputParams = zeros(1,numel(parameters));
    I_notFit = find(parameters == 0);
    for v = 1:numel(I_notFit)
        if I_notFit(v) == 1
            inputParams(1) = sigma_ellipse; % Use calibration data for default sigma_criterion
        elseif I_notFit(v) == 2
            inputParams(2) = 5; % Default sigma_criterion
        elseif I_notFit(v) == 3
            inputParams(3) = 0; % Default lapse (i.e., no lapse)
        elseif I_notFit(v) == 4
            inputParams(4) = Inf; % Default gamma
        elseif I_notFit(v) == 5
            inputParams(5) = .2; % Default alpha
        elseif I_notFit(v) == 6
            inputParams(6) = 1; % Default w (i.e., no bias)
        end
    end

    %% Choose priors - start with uniformative priors for all parameters
        % These will be uniform for all priors since we took the log of the
        % sigma parameters
    prior = ones(1, gridSize(1));
    multParams = ones(1, gridSize(1));
    for vv = 1:NumParams
        prior = unifpdf(params2fit(vv,:), min(params2fit(vv,:)), max(params2fit(vv,:)));
        prior = prior.*prior;
        multParams = multParams.*params2fit(vv,:);
    end
    minParams = min(multParams);
    maxParams = max(multParams);
    prior = prior/(qtrapz(prior)*(maxParams-minParams)/(numel(prior)-1));
    prior = log(prior);

    %% Compute the negative log likelihood for all parameter sets

    switch NumParams
        case 1
            nLL_mat = zeros(1,gridSize(1));
            rmse_mat = nLL_mat;
            resp_model_mat = zeros(NumTrials, gridSize(1));
            p_estimate_mat = resp_model_mat;
            for ii = 1:gridSize(1)
                inputParams(I_params) = exp(Grid(ii));
                [nLL_mat(ii), rmse_mat(ii), resp_model_mat(:,ii), p_estimate_mat(:,ii)] = changeprob_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, task, score, model);
            end
        case 2
            nLL_mat = zeros(gridSize(1), gridSize(2));
            rmse_mat = nLL_mat;
            resp_model_mat = cell(gridSize(1), gridSize(2));
            p_estimate_mat = resp_model_mat;
            for ii = 1:gridSize(1)
                for jj = 1:gridSize(2)
                    inputParams(I_params) = Grid{ii,jj};
                    if sum(I_params == 1) == 1
                        inputParams(1) = exp(inputParams(1));
                    elseif sum(I_params == 2) == 1
                        inputParams(2) = exp(inputParams(2));
                    end
                    [nLL_mat(ii,jj), rmse_mat(ii,jj), resp_model_mat{ii,jj}, p_estimate_mat{ii,jj}] = changeprob_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, task, score, model);
                end
            end
        case 3
            nLL_mat = zeros(gridSize(1), gridSize(2), gridSize(3));
            rmse_mat = nLL_mat;
            resp_model_mat = cell(gridSize(1), gridSize(2), gridSize(3));
            p_estimate_mat = resp_model_mat;
            for ii = 1:gridSize(1)
                for jj = 1:gridSize(2)
                    for kk = 1:gridSize(3)
                        inputParams(I_params) = Grid{ii,jj,kk};
                        if sum(I_params == 1) == 1
                            inputParams(1) = exp(inputParams(1));
                        elseif sum(I_params == 2) == 1
                            inputParams(2) = exp(inputParams(2));
                        end
                        [nLL_mat(ii,jj,kk), rmse_mat(ii,jj,kk), resp_model_mat{ii,jj,kk}, p_estimate_mat{ii,jj,kk}] = changeprob_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, task, score, model);
                    end
                end
            end 
        case 4
            nLL_mat = zeros(gridSize(1), gridSize(2), gridSize(3), gridSize(4));
            rmse_mat = nLL_mat;
            resp_model_mat = cell(gridSize(1), gridSize(2), gridSize(3), gridSize(4));
            p_estimate_mat = resp_model_mat;
            for ii = 1:gridSize(1)
                for jj = 1:gridSize(2)
                    for kk = 1:gridSize(3)
                        for qq = 1:gridSize(4)
                            inputParams(I_params) = Grid{ii,jj,kk,qq};
                            if sum(I_params == 1) == 1
                                inputParams(1) = exp(inputParams(1));
                            elseif sum(I_params == 2) == 1
                                inputParams(2) = exp(inputParams(2));
                            end
                            [nLL_mat(ii,jj,kk,qq), rmse_mat(ii,jj,kk,qq), resp_model_mat{ii,jj,kk,qq}, p_estimate_mat{ii,jj,kk,qq}] = changeprob_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, task, score, model);
                        end
                    end
                end
            end
        case 5
            nLL_mat = zeros(gridSize(1), gridSize(2), gridSize(3), gridSize(4), gridSize(5));
            rmse_mat = nLL_mat;
            resp_model_mat = cell(gridSize(1), gridSize(2), gridSize(3), gridSize(4), gridSize(5));
            p_estimate_mat = resp_model_mat;
            for ii = 1:gridSize(1)
                for jj = 1:gridSize(2)
                    for kk = 1:gridSize(3)
                        for qq = 1:gridSize(4)
                            for pp = 1:gridSize(5)
                                inputParams(I_params) = Grid{ii,jj,kk,qq,pp};
                                if sum(I_params == 1) == 1
                                    inputParams(1) = exp(inputParams(1));
                                elseif sum(I_params == 2) == 1
                                    inputParams(2) = exp(inputParams(2));
                                end
                                [nLL_mat(ii,jj,kk,qq,pp), rmse_mat(ii,jj,kk,qq,pp), resp_model_mat{ii,jj,kk,qq,pp}, p_estimate_mat{ii,jj,kk,qq,pp}] = changeprob_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, task, score, model);
                            end
                        end
                    end
                end
            end
    end

    %% Compute the marginal likelihood
    switch NumParams
        case 1
            % Add the log likelihood to the log prior, subtract the max, and
            % exponentiate
            modelPost = -nLL_mat + prior;
            modelPost = exp(modelPost-max(modelPost(:)));
            
            % Intergrate across model parameter
            %increment = params2fit(2)-params2fit(1);
            %marginalLikelihood = qtrapz(modelPost)*increment;
            marginalLikelihood = qtrapz(modelPost);

            % Find the best fitting parameter (MAP)
            [~, I] = max(modelPost);
            bestFit_param = params2fit(I);
            fitParams = exp(bestFit_param);

            % nLL, rmse, p_estimate, beta_model for best fit parameter
            nLL = nLL_mat(I);
            rmse = rmse_mat(I);
            p_estimate = p_estimate_mat(:,I);
            resp_model = resp_model_mat(:,I);
        case 2
            % Add the log likelihood to the log prior, subtract the max, and
            % exponentiate
            prior = repmat(prior, numel(prior), 1);
            modelPost = -nLL_mat + prior;
            modelPost = exp(modelPost - max(modelPost(:)));

            % Integrate across parameters
            %increment = params2fit(:,2)-params2fit(:,1);
            %marginalLikelihood = qtrapz(qtrapz(modelPost)*increment(1))*increment(2);
            marginalLikelihood = qtrapz(qtrapz(modelPost));

            % Find best fitting parameters (MAP)
            [~, I] = max(modelPost(:));
            [I_firstParam, I_secondParam] = ind2sub(size(modelPost), I);
            bestFit_firstParam = firstParam(I_firstParam);
            bestFit_secondParam = secondParam(I_secondParam);
            fitParams = [exp(bestFit_firstParam), bestFit_secondParam];

            % nLL, rmse, p_estimate, beta_model for best fit parameters
            nLL = nLL_mat(I_firstParam, I_secondParam);
            rmse = rmse_mat(I_firstParam, I_secondParam);
            p_estimate = p_estimate_mat{I_firstParam, I_secondParam};
            resp_model = resp_model_mat{I_firstParam, I_secondParam};
        case 3
            % Add the log likelihood to the log prior, subtract the max, and
            % exponentiate
            s_prior = numel(prior);
            prior = repmat(prior, [s_prior, 1, s_prior]);
            modelPost = -nLL_mat + prior;
            modelPost = exp(modelPost - max(modelPost(:)));

            % Integrate across parameters
            %increment = params2fit(:,2)-params2fit(:,1);
            %marginalLikelihood = qtrapz(qtrapz(qtrapz(modelPost)*increment(1))*increment(2))*increment(3);
            marginalLikelihood = qtrapz(qtrapz(qtrapz(modelPost)));

            % Find best fitting parameters (MAP)
            [~, I] = max(modelPost(:));
            [I_firstParam, I_secondParam, I_thirdParam] = ind2sub(size(modelPost), I);

            bestFit_firstParam = exp(firstParam(I_firstParam));
            bestFit_secondParam = secondParam(I_secondParam);
            bestFit_thirdParam = thirdParam(I_thirdParam);
            fitParams = [bestFit_firstParam, bestFit_secondParam, bestFit_thirdParam];

            % nLL, rmse, p_estimate, beta_model for best fit parameters
            nLL = nLL_mat(I_firstParam, I_secondParam, I_thirdParam);
            rmse = rmse_mat(I_firstParam, I_secondParam, I_thirdParam);
            p_estimate = p_estimate_mat{I_firstParam, I_secondParam, I_thirdParam};
            resp_model = resp_model_mat{I_firstParam, I_secondParam, I_thirdParam};
        case 4
            % Add the log likelihood to the log prior, subtract the max, and
            % exponentiate
            prior = repmat(prior, [numel(prior), 1, numel(prior), numel(prior)]);
            modelPost = -nLL_mat + prior;
            modelPost = exp(modelPost - max(modelPost(:)));

            % Integrate across parameters
            %increment = params2fit(:,2)-params2fit(:,1);
            %marginalLikelihood = qtrapz(qtrapz(qtrapz(qtrapz(modelPost)*increment(1))*increment(2))*increment(3))*increment(4);
            marginalLikelihood = qtrapz(qtrapz(qtrapz(qtrapz(modelPost))));

            % Find best fitting parameters (MAP)
            [~, I] = max(modelPost(:));
            [I_firstParam, I_secondParam, I_thirdParam, I_fourthParam] = ind2sub(size(modelPost), I);

            bestFit_firstParam = exp(firstParam(I_firstParam));
            bestFit_secondParam = secondParam(I_secondParam);
            bestFit_thirdParam = thirdParam(I_thirdParam);
            bestFit_fourthParam = fourthParam(I_fourthParam);
            fitParams = [bestFit_firstParam, bestFit_secondParam, bestFit_thirdParam, bestFit_fourthParam];

            % nLL, rmse, p_estimate, beta_model for best fit parameters
            nLL = nLL_mat(I_firstParam, I_secondParam, I_thirdParam, I_fourthParam);
            rmse = rmse_mat(I_firstParam, I_secondParam, I_thirdParam, I_fourthParam);
            p_estimate = p_estimate_mat{I_firstParam, I_secondParam, I_thirdParam, I_fourthParam};
            resp_model = resp_model_mat{I_firstParam, I_secondParam, I_thirdParam, I_fourthParam};
        case 5
            % Add the log likelihood to the log prior, subtract the max, and
            % exponentiate
            prior = repmat(prior, [numel(prior), 1, numel(prior), 1, numel(prior)]);
            modelPost = -nLL_mat + prior;
            modelPost = exp(modelPost - max(modelPost(:)));

            % Integrate across parameters
            %increment = params2fit(:,2)-params2fit(:,1);
            %marginalLikelihood = qtrapz(qtrapz(qtrapz(qtrapz(qtrapz(modelPost)*increment(1))*increment(2))*increment(3))*increment(4))*increment(5);
            marginalLikelihood = qtrapz(qtrapz(qtrapz(qtrapz(qtrapz(modelPost)))));

            % Find best fitting parameters (MAP)
            [~, I] = max(modelPost(:));
            [I_firstParam, I_secondParam, I_thirdParam, I_fourthParam, I_fifthParam] = ind2sub(size(modelPost), I);

            bestFit_firstParam = exp(firstParam(I_firstParam));
            bestFit_secondParam = secondParam(I_secondParam);
            bestFit_thirdParam = thirdParam(I_thirdParam);
            bestFit_fourthParam = fourthParam(I_fourthParam);
            bestFit_fifthParam = fifthParam(I_fifthParam);
            fitParams = [bestFit_firstParam, bestFit_secondParam, bestFit_thirdParam, bestFit_fourthParam, bestFit_fifthParam];

            % nLL, rmse, p_estimate, beta_model for best fit parameters
            nLL = nLL_mat(I_firstParam, I_secondParam, I_thirdParam, I_fourthParam, I_fifthParam);
            rmse = rmse_mat(I_firstParam, I_secondParam, I_thirdParam, I_fourthParam, I_fifthParam);
            p_estimate = p_estimate_mat{I_firstParam, I_secondParam, I_thirdParam, I_fourthParam, I_fifthParam};
            resp_model = resp_model_mat{I_firstParam, I_secondParam, I_thirdParam, I_fourthParam, I_fifthParam};
        otherwise
            error('changeprob_mL can only fit 1-5 parameters.')
    end
    toc
end