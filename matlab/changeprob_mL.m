function [marginalLikelihood, modelPost, nLL, rmse, fitParams, p_estimate,...
    p_true, resp_model, beta_resp] = changeprob_mL(model, data, task)

%% CHANGEPROB_ML Computes marginal likelihood for the specified model

    % Here we compute the marginal likelihood of the specfied model. 
    % To calculate the marginal likelihood, we multiply the 
    % likelihood of the data given the model by the prior of the model and 
    % integrate (i.e., sum) over the model parameters.
    
    % INPUT:
        % model: a character string indicating the model you want to fit, 
        % which is limited to the following:
            % 'idealBayesian'
            % 'exponential'
            % 'RL_probability'
            % 'exponential_conservative'
            % 'RL_probability_conservative'
            % 'RL_criterion'
            % 'fixed'
        % data: data struct from changing probability experiment
        % task: lets you choose which task to fit
            % 1 - overt-criterion task
            % 2 - covert-criterion task
        % options: work in progress
        
    % OUTPUT:
        % marginalLikelihood: a measure of model fit
        % modelPost: proportional to the p(parameters | data, model)
        % nLL: negative log likelihood that correponds to the best fitting
        % parameters
        % rmse: root mean squared error between the model probability and
        % the true probability (Note: does not apply to the RL_criterion model)
        % fitParams: best fitting model parameters computed by taking the
        % MAP of modelPost
        % p_estimate: model's estimate of probability (Note: does not apply to
        % the RL_criterion model)
        % p_true: true probability values
        % beta_model: predicted criterion
        % beta_resp: observer's criterion
        
    % Code adapted from Luigi Acerbi's changeprob_bocpd_nll.m 
    
    % Author:   Elyse norton
    % Email:    elyse.norton@gmail.com
    % Date:     10/4/2016
    
    % Model to be fit
    if nargin < 1; error('Please indicate the model you want to fit.'); end
    potentialModels = {'idealBayesian', 'exponential', 'RL_probability', ...
        'exponential_conservative', 'RL_probability_conservative', ...
        'RL_criterion', 'fixed'};
    model = find(strcmp(model, potentialModels)); % recode model to numeric value
    
    % Data struct or random seed for fake data generation
    if nargin < 2 || isempty(data); data = 0; end
    if isnumeric(data); rng(data); data = []; end

    % Task (1 overt, 2 covert)
    if nargin < 3 || isempty(task); task = 1; end % overt-criterion task is the default
    if task ~= 1 && task ~= 2; error('TASK can only be 1 (overt) or 2 (covert).'); end

    % Additional options to set experimental constants
    if nargin < 4; options = []; end

    %do_plot = nargout == 0; % If no outputs, make a plot
    do_plot = 1; % Always plot
    
    %% Read session parameters or create a fake dataset
    if ~isempty(data)
        NumTrials = data.NumTrials;

        col = data.SessionOrder(1);     % Column of overt criterion task
        
        sigma_ellipse = data.EllipseNoise;
        sigma_s = data.StdDev;
        sigma = sqrt(sigma_s^2 + sigma_ellipse^2);  % Category + sensory noise
        mu = [data.GreenMean(col),data.RedMean(col)];
        mu_bar = mean(mu);
        mu = mu - mu_bar;
        C = (data.Category(:,col) == 2); % Category A
        p_true = data.pA(:,col);
        catMeans = [data.RedMean(col),data.GreenMean(col)];
        S = data.StimulusAngle(:,col) - catMeans(data.Category(:,col))';
        % For the RL_criterion model we need to use noisy measurements of S
        X = zeros(NumTrials, 1000);
        for q = 1:1000
            X(:,q) = S + sigma_ellipse*randn(NumTrials, 1);
        end
        beta_resp = data.Criterion(:,col) - mu_bar;   % Reported criterion
        score = data.Score(:,col);
    else
        NumTrials = 800;
        sigma_s = 10;
        if isempty(sigma_ellipse) || ~isfinite(sigma_ellipse); sigma_ellipse = 5.68; end
        if isempty(sigma_criterion) || ~isfinite(sigma_criterion); sigma_criterion = sigma_ellipse; end
        sigma = sqrt(sigma_s^2 + sigma_ellipse^2);
        mu = [-10,10];

        C = []; p_true = []; p = 0;
        while numel(C) < NumTrials
            runlength = randi(prior_rl);
            pold = p;
            while p == pold; p = p_vec(randi(Nprobs)); end
            p_true = [p_true; p*ones(runlength,1)];
            C = [C; 2 - (rand(runlength,1) < p)];          % Vector of true categories
        end
        C = C(1:NumTrials);
        p_true = p_true(1:NumTrials);

        % Generate stimuli based on category labels
        S = mu(C)' + sigma_s*randn(NumTrials,1);

        % Generate noisy measurements for RL_criterion model
        X = zeros(NumTrials, 1000);
        for q = 1:size(X,2)
            X(:,q) = S + sigma_ellipse*randn(NumTrials, 1);
        end

        beta_resp = zeros(NumTrials,1);
    end
    
    %% Observer model parameters
    switch task
        case 1
            sigma_criterion = log(linspace(1,30));
            sigma_noise = sigma_criterion;
        case 2
            sigma_ellipse = log(linspace(1,30));
            sigma_noise = sigma_ellipse;
    end
    
    if model == 2 || model == 3 || model == 6
        eta = linspace(0,1); % 0 <= eta <= 1
        % Set up parameter grid
        s_eta = numel(eta);
        s_sigma_noise = numel(sigma_noise);
        Grid = cell(s_eta, s_sigma_noise);
        for i = 1:s_eta
            for j = 1:s_sigma_noise
                Grid{i,j} = [eta(i), sigma_noise(j)];
            end
        end
    elseif model == 4 || model == 5
        eta = linspace(0,1);
        w = linspace(0,1);
        % Set up parameter grid
        s_eta = numel(eta);
        s_sigma_noise = numel(sigma_noise);
        s_w = numel(w);
        Grid = cell(s_eta, s_sigma_noise, s_w);
        for i = 1:s_eta
            for j = 1:s_sigma_noise
                for k = 1:s_w
                    Grid{i,j,k} = [eta(i), sigma_noise(j), w(k)];
                end
            end
        end
    end

    %% Choose priors
    prior_sigma_noise = unifpdf(sigma_noise, min(sigma_noise), max(sigma_noise)); % Jefferys prior - uniform in log space
    
    if model == 2 || model == 3 || model == 6
        % Start with uniformative priors
        prior_eta = unifpdf(eta, 0, 1); % Uniform prior

        % Multiply and scale
        prior = prior_eta.*prior_sigma_noise;
        prior = prior/(trapz(prior)*(max(eta.*sigma_noise)-min(eta.*sigma_noise))/(numel(prior)-1));
        prior = log(prior);
    elseif model == 4 || model == 5
        % Start with uniformative priors
        prior_eta = unifpdf(eta, 0, 1); % Uniform prior
        prior_w = unifpdf(w, 0, 1); % Uniform prior

        % Multiply and scale
        prior = prior_eta.*prior_sigma_noise.*prior_w;
        prior = prior/(trapz(prior)*(max(eta.*sigma_noise.*w)-min(eta.*sigma_noise.*w))/(numel(prior)-1));
        prior = log(prior);
    elseif model == 1 || model == 7
        prior = log(prior_sigma_noise);
    end

    %% Compute the negative log likelihood for all parameter sets

    if model == 2 || model == 3 || model == 6
        nLL_mat = zeros(numel(eta), numel(sigma_noise));
        rmse_mat = zeros(numel(eta), numel(sigma_noise));
        resp_model_mat = cell(size(Grid));
        p_estimate_mat = cell(size(Grid));
        for ii = 1:numel(eta)
            for jj = 1:numel(sigma_noise)
                parameters = Grid{ii,jj};
                parameters(2) = exp(parameters(2));
                if model == 2 || model == 3
                    [nLL_mat(ii,jj), resp_model_mat{ii,jj}, rmse_mat(ii,jj), p_estimate_mat{ii,jj}] = changeprob_nLL(parameters);
                elseif model == 6
                    [nLL_mat(ii,jj), resp_model_mat{ii,jj}] = changeprob_nLL(parameters);
                end
            end
        end
    elseif model == 4 || model == 5
        nLL_mat = zeros(numel(eta), numel(sigma_noise), numel(w));
        rmse_mat = zeros(numel(eta), numel(sigma_noise), numel(w));
        resp_model_mat = cell(size(Grid));
        p_estimate_mat = cell(size(Grid));
        for ii = 1:numel(eta)
            for jj = 1:numel(sigma_noise)
                for kk = 1:numel(w) 
                    parameters = Grid{ii,jj,kk};
                    parameters(2) = exp(parameters(2));
                    [nLL_mat(ii,jj,kk), resp_model_mat{ii,jj,kk}, rmse_mat(ii,jj,kk), p_estimate_mat{ii,jj,kk}] = changeprob_nLL(parameters);
                end
            end
        end
    elseif model == 1 || model == 7
        nLL_mat = zeros(numel(sigma_noise), 1);
        rmse_mat = zeros(numel(sigma_noise), 1);
        resp_model_mat = zeros(numel(sigma_noise), NumTrials);
        p_estimate_mat = zeros(numel(sigma_noise), NumTrials);
        for ii = 1:numel(sigma_noise)
            parameters = exp(sigma_noise(ii));
            [nLL_mat(ii), resp_model_mat(ii,:), rmse_mat(ii), p_estimate_mat(ii,:)] = changeprob_nLL(parameters);
        end
    end

    %% Compute the model posterior
    switch numel(parameters)
        case 1
            % Add the log likelihood to the log prior, subtract the max, and
            % exponentiate
            modelPost = -nLL_mat' + prior;
            modelPost = exp(modelPost-max(modelPost(:)));

            % Intergrate across model parameter
            marginalLikelihood = trapz(modelPost);

            % Find the best fitting parameter (MAP)
            [~, I] = max(modelPost);
            bestFit_SigmaNoise = sigma_noise(I);
            fitParams = exp(bestFit_SigmaNoise);

            % nLL, rmse, p_estimate, beta_model for best fit parameter
            nLL = nLL_mat(I);
            rmse = rmse_mat(I);
            p_estimate = p_estimate_mat(I,:);
            resp_model = resp_model_mat(I,:);
        case 2
            % Add the log likelihood to the log prior, subtract the max, and
            % exponentiate
            prior = repmat(prior, numel(prior), 1);
            modelPost = -nLL_mat + prior;
            modelPost = exp(modelPost - max(modelPost(:)));

            % Integrate across parameters
            marginalLikelihood = trapz(trapz(modelPost));

            % Find best fitting parameters (MAP)
            [~, I] = max(modelPost(:));
            [I_eta, I_sigma_noise] = ind2sub(size(modelPost), I);
            bestFit_Eta = eta(I_eta);
            bestFit_SigmaNoise = sigma_noise(I_sigma_noise);
            fitParams = [bestFit_Eta, exp(bestFit_SigmaNoise)];

            % nLL, rmse, p_estimate, beta_model for best fit parameters
            nLL = nLL_mat(I_eta, I_sigma_noise);
            rmse = rmse_mat(I_eta, I_sigma_noise);
            p_estimate = p_estimate_mat{I_eta, I_sigma_noise};
            resp_model = resp_model_mat{I_eta, I_sigma_noise};
        case 3
            % Add the log likelihood to the log prior, subtract the max, and
            % exponentiate
            s_prior = numel(prior);
            prior = repmat(prior, [s_prior, 1, s_prior]);
            modelPost = -nLL_mat + prior;
            modelPost = exp(modelPost - max(modelPost(:)));

            % Integrate across parameters
            marginalLikelihood = trapz(trapz(trapz(modelPost)));

            % Find best fitting parameters (MAP)
            [~, I] = max(modelPost(:));
            [I_eta, I_sigma_noise, I_w] = ind2sub(size(modelPost), I);

            bestFit_Eta = eta(I_eta);
            bestFit_SigmaNoise = sigma_noise(I_sigma_noise);
            bestFit_w = w(I_w);
            fitParams = [bestFit_Eta, exp(bestFit_SigmaNoise), bestFit_w];

            % nLL, rmse, p_estimate, beta_model for best fit parameters
            nLL = nLL_mat(I_eta, I_sigma_noise, I_w);
            rmse = rmse_mat(I_eta, I_sigma_noise, I_w);
            p_estimate = p_estimate_mat{I_eta, I_sigma_noise, I_w};
            resp_model = resp_model_mat{I_eta, I_sigma_noise, I_w};
    end

%% Plots
    if do_plot
        
        if model ~= 4 && model ~= 5
            % Model posterior - 1 parameter
            if model == 1 || model == 7
                figure(1); hold on;
                plot(sigma_noise, modelPost, '-k');
                xlabel('sigma_noise');
                ylabel('Posterior density');
                hold off;
            else
                % Model posterior (surface) - 2 parameter
                figure(1); hold on;
                surf(sigma_noise, eta, modelPost);
                xlabel('sigma_noise');
                ylabel('Eta');
                zlabel('Posterior density');
                hold off;
            end
        end
                
        switch task
            case 1 
                if model ~= 6 && model ~= 1 
                    % p_estimate vs. p_true and resp_model vs. beta_resp
                    figure(2);
                    subplot(2,1,1); hold on;
                    plot(1:NumTrials, p_true, '-k', 'LineWidth', 2);
                    plot(1:NumTrials, smooth(p_estimate,11), '-b', 'LineWidth', 2);
                    ylim([0,1]);
                    box off;
                    set(gca,'XTickLabel',[]);
                    ylabel('Pr(A)');    
                    text(NumTrials*0.01,0.9,['RMSE = ' num2str(rmse,'%.3g')]);

                    % Smooth observed criterion for visualization
                    subplot(2,1,2);
                    beta_smooth = smooth(beta_resp,11);
                    scatter(1:NumTrials,beta_resp,'k.'); hold on;
                    plot(1:NumTrials,beta_smooth,'k','LineWidth',2); hold on;
                    plot(1:NumTrials,smooth(resp_model,11),'r','LineWidth',2);
                    ylabel('Criterion');
                    ylim([-50 120]);
                    set(gca,'YTick',[-50,0,50]);
                    h_leg = legend('Reported criterion','Moving average','Model');
                    set(h_leg,'Box','off','Location','NorthEast');
                    xlabel('Trial number');
                    hold off;
                elseif model == 6
                    figure(2); hold on;
                    beta_smooth = smooth(beta_resp,11);
                    scatter(1:NumTrials,beta_resp,'k.'); hold on;
                    plot(1:NumTrials,beta_smooth,'k','LineWidth',2); hold on;
                    plot(1:NumTrials,smooth(resp_model,11),'r','LineWidth',2);
                    ylabel('Criterion');
                    ylim([-50 120]);
                    set(gca,'YTick',[-50,0,50]);
                    h_leg = legend('Reported criterion','Moving average','Model');
                    set(h_leg,'Box','off','Location','NorthEast');
                    xlabel('Trial number');
                    hold off;
%                 else
%                     h = plotify([1 1, 2:Nprobs+1, (Nprobs+2)*[1 1 1]]','Margins',[0.05 0.05, 0.15 0.05],'Title',['Posterior over category probability for the BOCPD ideal observer']);
%                     axes(h(1));
%                     plot(1:NumTrials,p_true,'k','LineWidth',2); hold on;
%                     plot(1:NumTrials,meanP,'r','LineWidth',2);
%                     ylim([0,1]);
%                     box off;
%                     set(gca,'XTickLabel',[]);
%                     ylabel('Pr(A)');    
%                     text(NumTrials*0.01,0.9,['RMSE = ' num2str(rmse,'%.3g')]);
% 
%                     change_vec = [0;find(diff(p_true) ~= 0);NumTrials];
% 
%                     for i = 1:Nprobs
%                         axes(h(i+1));
% 
%                         for j = 1:numel(change_vec)-1
%                             plot(change_vec(j)*[1 1],[0 1],':k','LineWidth',1);
%                             if abs(p_true(change_vec(j)+1) - p_vec(i)) < 1e-8
%                             % if p_true(change_vec(j)+1) == p_vec(i)
%                                 xx = [change_vec(j) change_vec(j+1) change_vec(j+1) change_vec(j)];
%                                 yy = [0 0, 1 1];
%                                 % fill(xx,yy,'c','EdgeColor','none','FaceAlpha',0.5);                
%                                 fill(xx,yy,'c','EdgeColor','none'); hold on;
%                             end
%                         end
% 
%                         plot(1:NumTrials,P(1:NumTrials,i),'k'); hold on;
%                         box off;
%                         if i < Nprobs; set(gca,'XTickLabel',[]); end
%                         text(NumTrials*0.01,1.1,['Pr(A) = ' num2str(p_vec(i))])
%                     end
% 
%                     axes(h(Nprobs+2));
% 
%                     % Smoothen responded criterion for visualization
%                     beta_smooth = smooth(beta_resp,11);
% 
%                     if ~isempty(data)
%                         scatter(1:NumTrials,beta_resp,'k.'); hold on;
%                         plot(1:NumTrials,beta_smooth,'k','LineWidth',2); hold on;
%                     end
%                     plot(1:NumTrials,beta_opt(1:NumTrials),'r','LineWidth',2);
%                     ylabel('Criterion');
%                     ylim([-50 120]);
%                     set(gca,'YTick',[-50,0,50]);
%                     if ~isempty(data)
%                         h_leg = legend('Reported criterion','Moving average','Model');
%                     else
%                         h_leg = legend('Model');                
%                     end
%                     set(h_leg,'Box','off','Location','NorthEast');
%                     xlabel('Trial number');
                    error('Ideal Bayesian model not yet fully incorporated');
                end
            case 2
                error('Covert-criterion task not coded yet!');
        end
    end

    %-----------------------------------------------------------------------------------------
    function [nLL, model_resp, rmse, pi_estimate] = changeprob_nLL(parameters, p_initial, p_bias)
    %CHANGEPROB_NLL computes the negative log likelihood for data in the changing
    %probabilities experiment given a specified model

        % ALL model parameters must be specified
        if nargin < 1; error('You must specify parameter values.'); end

        if and(model == 2, numel(parameters) ~= 2) || and(model == 3, numel(parameters) ~= 2)...
                || and(model == 6, numel(parameters) ~= 2)
            errorText = strcat('The', {' '}, potentialModels{model}, {' '}, 'model must have exactly 2 parameters.');
            error(errorText);
        elseif and(model == 1, numel(parameters ~= 1))
            error('The ideal Bayesian model must have exactly 1 parameter.');
        elseif and(model == 4, numel(parameters) ~= 3) || and(model == 5, numel(parameters) ~= 3)
            errorText = strcat('The', {' '}, potentialModels{model}, {' '}, 'model must have exactly 3 parameters.');
            error(errorText);
        end

        % Initial probability estimate
        if nargin < 2 || isempty(p_initial); p_initial = .5; end % default to equal probability

        % Initial probability estimate
        if nargin < 3 || isempty(p_bias); p_bias = .5; end % default to .5

        % INPUT

            % model: specify an observer model from the following list:
                % 'idealBayesian'
                % 'exponential'
                % 'RL_probability'
                % 'RL_criterion'
                % 'exponential_conservative'
                % 'RL_probability_conservative'
            % parameters: vector of free parameters
                  % For the ideal Bayesian model
                    % parameters(1): adjustment (overt) or sensory (covert) noise
                  % For the exponential, RL_probability, or RL_criterion
                  % models
                    % parameters(1): smoothing factor/learning rate
                    % parameters(2): adjustment (overt) or sensory (covert) noise
                  % For the exponential_conservative and the RL_probability_conservative
                    % parameters(1): smoothing factor/learning rate
                    % parameters(2): adjustment (overt) or sensory (covert) noise
                    % parameters(3): weight given to probability estimate
            % p_initial: starting probability assumed by observer (default is .5)
            % p_bias: probability bias

        % OUTPUT

            % nLL: negative log likelihood
            % model_resp: model criterion (overt) OR model categorizations (covert)
            % rmse: root mean squared error between model estimates of pi and p_true
                % NOTE: Not applicable to the RL_criterion model
            % pi_estimate: vector of the model's estimates of pi
                % NOTE: Not applicable to the RL_criterion model

        n = numel(C);
        log_Pbeta = zeros(n,1);
        model_resp = zeros(n,1);
        if task == 1        
            if model == 2
                Gamma = zeros(n,1);
                pi_estimate = zeros(n,1);
                pi_estimate(1,:) = p_initial;
                for t = 1:n
                    Gamma(t) = pi_estimate(t,:)/(1-pi_estimate(t,:));
                    %model_resp(t) = sigma^2 * log(Gamma(t)) / diff(mu) + parameters(2)*randn(1,1);
                    model_resp(t) = sigma^2 * log(Gamma(t)) / diff(mu); % This in incorrect need to include noise in model
                    log_Pbeta(t) = -0.5*log(2*pi*parameters(2)) - 0.5*((beta_resp(t)-model_resp(t))./parameters(2)).^2;
                    pi_estimate(t+1,:) = pi_estimate(t,:) + parameters(1)*(C(t)-pi_estimate(t,:));
                end

                pi_estimate = pi_estimate(1:numel(C));
                nLL = -nansum(log_Pbeta);
                rmse = sqrt(mean((pi_estimate - p_true).^2));
            elseif model == 3
                Gamma = zeros(n,1);
                pi_estimate = zeros(n,1);
                pi_estimate(1,:) = p_initial;
                for t = 1:n
                    Gamma(t) = pi_estimate(t,:)/(1-pi_estimate(t,:));
                    %model_resp(t) = sigma^2 * log(Gamma(t)) / diff(mu) + parameters(2)*randn(1,1);
                    model_resp(t) = sigma^2 * log(Gamma(t)) / diff(mu);
                    log_Pbeta(t) = -0.5*log(2*pi*parameters(2)) - 0.5*((beta_resp(t)-model_resp(t))./parameters(2)).^2;
                    if score(t) == 0
                        pi_estimate(t+1,:) = pi_estimate(t,:) + parameters(1)*(C(t)-pi_estimate(t,:));
                    else
                        pi_estimate(t+1,:) = pi_estimate(t,:);
                    end
                end

                pi_estimate = pi_estimate(1:numel(C));
                nLL = -nansum(log_Pbeta);
                rmse = sqrt(mean((pi_estimate - p_true).^2));
            elseif model == 4
                Gamma = zeros(n,1);
                pi_conservative = zeros(n,1);
                pi_estimate(1,:) = p_initial;
                for t = 1:n
                    pi_conservative(t,:) = parameters(3)*pi_estimate(t,:) + (1-parameters(3))*p_bias;
                    Gamma(t) = pi_conservative(t,:)/(1-pi_conservative(t,:));
                    %model_resp(t) = sigma^2 * log(Gamma(t)) / diff(mu) + parameters(2)*randn(1,1);
                    model_resp(t) = sigma^2 * log(Gamma(t)) / diff(mu);
                    log_Pbeta(t) = -0.5*log(2*pi*parameters(2)) - 0.5*((beta_resp(t)-model_resp(t))./parameters(2)).^2;
                    pi_estimate(t+1,:) = pi_estimate(t,:) + parameters(1)*(C(t)-pi_estimate(t,:));
                end

                pi_estimate = pi_conservative;
                nLL = -nansum(log_Pbeta);
                rmse = sqrt(mean((pi_conservative - p_true).^2));
            elseif model == 5
                Gamma = zeros(n,1);
                pi_conservative = zeros(n,1);
                pi_estimate(1,:) = p_initial;
                for t = 1:n
                    pi_conservative(t,:) = parameters(3)*pi_estimate(t,:) + (1-parameters(3))*p_bias;
                    Gamma(t) = pi_conservative(t,:)/(1-pi_conservative(t,:));
                    %model_resp(t) = sigma^2 * log(Gamma(t)) / diff(mu) + parameters(2)*randn(1,1);
                    model_resp(t) = sigma^2 * log(Gamma(t)) / diff(mu);
                    log_Pbeta(t) = -0.5*log(2*pi*parameters(2)) - 0.5*((beta_resp(t)-model_resp(t))./parameters(2)).^2;
                    if score(t) == 0
                        pi_estimate(t+1,:) = pi_estimate(t,:) + parameters(1)*(C(t)-pi_estimate(t,:));
                    else
                        pi_estimate(t+1,:) = pi_estimate(t,:);
                    end
                end

                pi_estimate = pi_conservative;
                nLL = -nansum(log_Pbeta);
                rmse = sqrt(mean((pi_conservative - p_true).^2));
            elseif model == 6
                % For each vector of noisy measurements fit model and average over all vectors
                nLL = zeros(1,size(X,2));
                for qq = 1:size(X,2)
                    x = X(:,qq); % vector of noisy measurements
                    Gamma(1) = p_initial/(1-p_initial);
                    %model_resp(1,qq) = sigma^2 * log(Gamma(1)) / diff(mu) + parameters(2)*randn(1,1);
                    model_resp(1,qq) = sigma^2 * log(Gamma(1)) / diff(mu);
                    for t = 1:n
                        log_Pbeta(t,qq) = -0.5*log(2*pi*parameters(2)) - 0.5*((beta_resp(t)-model_resp(t,qq))./parameters(2)).^2;
                        if score(t) == 0
                            %model_resp(t+1,qq) = model_resp(t,qq) + parameters(1)*(x(t)-model_resp(t,qq)) + parameters(2)*randn(1,1);
                            model_resp(t+1,qq) = model_resp(t,qq) + parameters(1)*(x(t)-model_resp(t,qq));
                        else
                            %model_resp(t+1,qq) = model_resp(t,qq) + parameters(2)*randn(1,1);
                            model_resp(t+1,qq) = model_resp(t,qq);
                        end
                    end
                    nLL(:,qq) = -nansum(log_Pbeta(:,qq));
                end
                model_resp = mean(model_resp,2);
                model_resp = model_resp(1:n);
                nLL = mean(nLL);
            elseif model == 7
                pi_estimate = zeros(n,1) + p_initial; % Assume probability is fixed throughout the experimental block
                Gamma = pi_estimate./(1-pi_estimate);
                %model_resp = sigma^2 * log(Gamma)/diff(mu) + parameters(1)*randn(n,1);
                model_resp = sigma^2 * log(Gamma)/diff(mu);
                log_Pbeta = -0.5*log(2*pi*parameters(1)) - 0.5*((beta_resp-model_resp)./parameters(1)).^2;
                
                nLL = -nansum(log_Pbeta);
                rmse = sqrt(mean((pi_estimate - p_true).^2));
            elseif model == 1
                error('Ideal Bayesian model not incorporated yet!');
            end
        else
            error('Covert-criterion task not available yet!');
        end
    end
end