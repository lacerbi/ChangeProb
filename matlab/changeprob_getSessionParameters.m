function [NumTrials, sigma_ellipseData, mu, sigma_s, C, S, p_true, resp, score] = changeprob_getSessionParameters(data, task, parameters)
%CHANGEPROB_GETSESSIONPARAMETERS Gets session parameters from an existing
%data struct or creates a fake dataset

%   INPUT: 
        % data
            % Experimental data struct to be decomposed
        % task: 1 - overt (default), 2 - covert, 3 - mixed
        % parameters used to generate fake data
            % parameters(1): sensory noise (sigma_ellipse)
            % parameters(2): adjustment noise (sigma_criterion)
        
%   OUTPUT:
        % NumTrials: total number of trials
        % sigma_ellipseData: sensory noise from calibration data
        % mu: vector containing the category means [muA, muB]
        % sigma_s: std dev of the category distributions
        % C: vector of category values (1 - A, 0 - B)
        % S: vector of true stimulus orientations
        % p_true: vector containing the probability of A
        % resp: vector containing the observer's criterion setting 
        % (overt task), categorizations (covert), or both (mixed)
        % score: 0 - wrong, 1 - correct
        
%   Author: Elyse Norton
%   Date: 1/11/17
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
        if task == 3
            NumTrials = data.NumTrials;    
            
            % Noise parameters
            sigma_s = data.StdDev;
            if isempty(sigma_ellipse) || ~isfinite(sigma_ellipse)
                sigma_ellipseData = data.EllipseNoise;
                sigma = sqrt(sigma_s^2 + sigma_ellipseData^2);
            else
                sigma = sqrt(sigma_s^2 + sigma_ellipse^2);
            end

            % Category information 
            % (in the data Category B is coded as 1, Category A is coded as 2)
            C = (data.TrialType == 2);        % 1 - A, 0 - B
            C = double(C);
            p_true = data.pA;
            mu = [data.MeanSignal,data.MeanNoise];

            % Shift coordinate system to zero
            mu_bar = mean(mu);  
            mu = bsxfun(@minus, mu, mu_bar);
            S = bsxfun(@minus, data.TrueAngle, mu_bar); 

            % Get task-relevant responses
            resp = data.response; 
            % Shift coordinate system for criterion responses
            I_overt = 5:5:NumTrials;
            resp(I_overt) = bsxfun(@minus, resp(I_overt), mu_bar);
            % Recode category response such that 1 - A, 0 - B
            resp_covert = resp;
            resp_covert(5:5:NumTrials) = [];
            resp_covert = double(resp_covert == 2);
            I_covert = 1:NumTrials;
            I_covert(5:5:NumTrials) = [];
            resp(I_covert) = resp_covert;
            score = data.score;
        else  
            col = data.SessionOrder(task);     % Column of overt-criterion task
            NumTrials = data.NumTrials;    

            % Noise parameters
            sigma_s = data.StdDev;
            if isempty(sigma_ellipse) || ~isfinite(sigma_ellipse)
                sigma_ellipseData = data.EllipseNoise; 
                sigma = sqrt(sigma_s^2 + sigma_ellipseData^2);
            else
                sigma = sqrt(sigma_s^2 + sigma_ellipse^2);
            end

            % Category information 
            % (in the data Category B/Red is coded as 1, Category A/Green is coded as 2)
            C = (data.Category(:,col) == 2);        % Category A/Green
            C = double(C);
            p_true = data.pA(:,col);
            mu = [data.GreenMean(col),data.RedMean(col)];

            % Shift coordinate system to zero
            mu_bar = mean(mu);
            mu = bsxfun(@minus, mu, mu_bar);
            S = bsxfun(@minus, data.StimulusAngle(:,col), mu_bar);

            % Get task-relevant responses
            switch task
                case 1  % Overt-criterion task    
                    resp = bsxfun(@minus, data.Criterion(:,col), mu_bar);   % Reported criterion
                case 2  % Covert-criterion task
                    resp = data.Response(:,col) == 2;
                    resp = double(resp);
            end
            score = data.Score(:,col);
        end   
    else  
        % Create fake data
        NumTrials = 800;
        sigma_s = 10;
        mu = [-8,8];
        prior_rl = [80,120];
        p_vec = linspace(0.2,0.8,5);
        Nprobs = numel(p_vec);
        if isempty(sigma_ellipse) || ~isfinite(sigma_ellipse)
            sigma_ellipseData = 5;
            sigma = sqrt(sigma_s^2 + sigma_ellipseData^2);
        else
            sigma = sqrt(sigma_s^2 + sigma_ellipse^2);
        end
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
        S = bsxfun(@plus, mu(C)', sigma_s*randn(NumTrials,1));
        if isempty(sigma_ellipse) || ~isfinite(sigma_ellipse)
            X = bsxfun(@plus, S, sigma_ellipseData*randn(NumTrials,1));
        else
            X = bsxfun(@plus, S, sigma_ellipse*randn(NumTrials,1));
        end

        % Responses based on fixed criterion (fixed at the neutral criterion)
        z_resp = bsxfun(@plus, mean(mu), sigma_criterion*randn(NumTrials,1));
        resp(:,1) = z_resp;

        Chat = ones(NumTrials,1);
        I = find(X > mean(mu));
        Chat(I) = Chat(I)+1;
        resp(:,2) = Chat;

        resp = resp(:,task);

        % Score
        scoreOvert = zeros(NumTrials, 1);
        I = find(or(and(C == 1, S <= z_resp), and(C == 2, S > z_resp)));
        scoreOvert(I) = bsxfun(@plus, scoreOvert(I), 1);
        score(:,1) = scoreOvert;

        scoreCovert = zeros(NumTrials, 1);
        I = find(C == Chat);
        scoreCovert(I) = scoreCovert(I)+1;
        score(:,2) = scoreCovert;

        score = score(:,task);
        C = double(C==1);
    end
end

