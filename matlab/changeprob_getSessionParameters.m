function [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp, score] = changeprob_getSessionParameters(data, task, parameters)
%CHANGEPROB_GETSESSIONPARAMETERS Gets session parameters from an existing
%data struct or creates a fake dataset

%   INPUT: 
        % data
            % Experimental data struct to be decomposed
        % task: 1 - overt (default), 2 - covert
        % parameters used to generate fake data
            % parameters(1): sensory noise (sigma_ellipse)
            % parameters(2): adjustment noise (sigma_criterion)
%   OUTPUT:
        % NumTrials: total number of trials
        % sigma_ellipse: sensory noise from calibration data
        % mu: vector containing the category means [muA, muB]
        % sigma: std dev of the internal distributions - sqrt(sigma_s^2 +
        % sigma_v^2)
        % C: vector of category values
        % S: vector of true stimulus angles
        % p_true: vector containing the probability of A
        % resp: vector containing the fake observer's criterion setting 
        % (overt task) or categorizations
        % score: 0 - wrong, 1 - correct
        
%   Author: Elyse Norton
%   Date: 10/18/16
%   email: elyse.norton@gmail.com
        
    switch nargin
        case 0
            data = [];
            task = 1;
            parameters = [];
        case 1
            task = 1;
            parameters = [];
        case 2
            parameters = [];
    end

    switch numel(parameters)
        case 0  
            sigma_ellipse = []; % Use calibration data
            sigma_criterion = []; % Use calibration data
        case 1
            sigma_ellipse = parameters(1);
        otherwise 
            sigma_ellipse = parameters(1);
            sigma_criterion = parameters(2);
    end
    
    if ~isempty(data)
        col = data.SessionOrder(task);     % Column of overt-criterion task
        NumTrials = data.NumTrials;    

        % Noise parameters
        sigma_s = data.StdDev;
        if isempty(sigma_ellipse) || ~isfinite(sigma_ellipse); sigma_ellipse = data.EllipseNoise; end
        sigma = sqrt(sigma_s^2 + sigma_ellipse^2);

        % Category information 
        % (in the data Category B/Red is coded as 1, Category A/Green is coded as 2)
        C = (data.Category(:,col) == 2);        % Category A/Green
        C = double(C);
        p_true = data.pA(:,col);
        mu = [data.GreenMean(col),data.RedMean(col)];

        % Shift coordinate system to zero
        mu_bar = mean(mu);
        mu = mu - mu_bar;
        S = data.StimulusAngle(:,col) - mu_bar;

        % Get task-relevant responses
        switch task
            case 1  % Overt-criterion task    
                resp = data.Criterion(:,col) - mu_bar;   % Reported criterion
            case 2  % Covert-criterion task
                resp = data.Response(:,col) == 2;
                resp = double(resp);
        end
        score = data.Score(:,col);
    else  
        % Create fake data
        NumTrials = 800;
        sigma_s = 10;
        mu = [-8,8];
        prior_rl = [80,120];
        p_vec = linspace(0.2,0.8,5);
        Nprobs = numel(p_vec);
        if isempty(sigma_ellipse) || ~isfinite(sigma_ellipse); sigma_ellipse = 5; end
        if isempty(sigma_criterion) || ~isfinite(sigma_criterion); sigma_criterion = sigma_ellipse; end
        sigma = sqrt(sigma_s^2 + sigma_ellipse^2);
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
        X = S + sigma_ellipse*randn(NumTrials,1);

        % Responses based on fixed criterion (fixed at the neutral criterion)
        z_resp = mean(mu) + sigma_criterion*randn(NumTrials,1);
        resp(:,1) = z_resp;

        Chat = ones(NumTrials,1);
        I = find(X > mean(mu));
        Chat(I) = Chat(I)+1;
        resp(:,2) = Chat;

        resp = resp(:,task);

        % Score
        scoreOvert = zeros(NumTrials, 1);
        I = find(or(and(C == 1, S <= z_resp), and(C == 2, S > z_resp)));
        scoreOvert(I) = scoreOvert(I)+1;
        score(:,1) = scoreOvert;

        scoreCovert = zeros(NumTrials, 1);
        I = find(C == Chat);
        scoreCovert(I) = scoreCovert(I)+1;
        score(:,2) = scoreCovert;

        score = score(:,task);
        C = double(C==1);
    end
end

