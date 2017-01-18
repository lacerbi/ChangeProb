function dataSim = changeprob_simulateData(subID, models, task, numSims)
%CHANGEPROB_SIMULATEDATA Simulates an observer with the same experimental
%parameters of the specified subject(s) using the specified model(s) in the
%specified task(s) numSims times.

% (Detailed documentation to be written.)
%
% Author:   Elyse Norton
% Email:    elyse.norton@gmail.com
% Date:     Jan/9/2017

if nargin < 1; subID = []; models = []; task = []; numSims = []; end

if isempty(subID)
    subID = {'CWG', 'EGC', 'EHN', 'ERK', 'GK', 'HHL', 'JKT', 'JYZ', 'RND', 'SML', 'SQC'}; % Default all subjects for tasks 1 and 2
    subID_mixed = {'CWG', 'EGC', 'EHN', 'ERK', 'HHL', 'RND', 'SML'}; % Default all subjects for task 3
end

if nargin < 2 || isempty(models)
    models = {'fixed', 'idealBayesian', 'exponential', 'RL_probability', ...
    'exponential_conservative', 'RL_probability_conservative', 'RL_criterion', ...
    'subBayesian_rlprior', 'subBayesian_conservative', 'subBayesian_pVec'}; % Default simulate all models
end

if nargin < 3 || isempty(task); task = [1 2 3]; end % Default simulate all tasks

if nargin < 4 || isempty(numSims); numSims = 1; end % Default number of simulations is 1/person/task & model

paramBounds_def = [1,30; 1,30; 0,0.1; -Inf,Inf; 0,1; 0,1; 2,200; 0,.5]; % Default parameter bounds
gridSize = 100; % Default grid size

%% For each subject, model, and task simulate an observer

currentFolder = pwd;

for i = 1:numel(task)
    if task(i) == 1
        currentTask = {'Overt'};
    elseif task(i) == 2
        currentTask = {'Covert'};
    else
        currentTask = {'Mixed'};
    end
    if task(i) ~= 3
        for j = 1:numel(models)
            if j < 3
                if i == 1
                    parameters = [0 1 0 0 0 0 0 0];
                else
                    parameters = [1 0 0 0 0 0 0 0];
                end
            elseif j == 3 || j == 4 || j == 7
                if i == 1
                    parameters = [0 1 0 0 1 0 0 0];
                else
                    parameters = [1 0 0 0 1 0 0 0];
                end
            elseif j == 8
                if i == 1
                    parameters = [0 1 0 0 0 0 1 0];
                else
                    parameters = [1 0 0 0 0 0 1 0];
                end
            elseif j == 9
                if i == 1
                    parameters = [0 1 0 0 0 1 0 0];
                else
                    parameters = [1 0 0 0 0 1 0 0];
                end
            elseif j == 10
                if i == 1
                    parameters = [0 1 0 0 0 0 0 1];
                else
                    parameters = [1 0 0 0 0 0 0 1];
                end
            else
                if i == 1
                    parameters = [0 1 0 0 1 1 0 0];
                else
                    parameters = [1 0 0 0 1 1 0 0];
                end
            end
            NumParams = sum(parameters);
            I_params = find(parameters ~= 0);
            paramBounds = paramBounds_def(I_params,:);
            for iParam = 1:NumParams
                switch I_params(iParam)
                    case {1, 2}
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize); % sigma_noise
                    case 3
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize); % lambda
                    case 4
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize); % gamma
                    case 5
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize); % alpha
                    case 6
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize); % w
                    case 7
                        params2fit(iParam,:) = round(linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize)); % Tmax (discrete)
                    case 8
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize); % Range for minimum probability (pVec(1))
                end
            end
            for ii = 1:numel(subID)
                % Load data
                cd('/Users/elysenorton/Desktop/ChangeProb/data');
                load(strcat('ChangingProbabilities_', subID{ii}));
                % Load fit data
                cd('/Users/elysenorton/Desktop/ChangeProb/matlab/ModelFit_data/SMALLlapse');
                load(char(strcat(subID{ii}, '_', models{j}, '_', currentTask)));
                cd(currentFolder);
                % Choose parameter vector with probability proportional to the joint posterior over parameters for given subject
                postVector = modelPost(:)/sum(modelPost(:));
                % Sort posterior vector
                [postSorted, I_post] = sortrows(postVector);
                % Compute the cumulative sum of sorted vector
                cumPost = cumsum(postSorted);
                for jj = 1:numSims
                    dataSim.randomSeed(jj,:) = i*j*ii*jj;
                    rng(dataSim.randomSeed(jj,:));
                    % Sample from posterior
                    I_cumsum = find(rand(1,1) < cumPost);
                    I_cumsum = I_cumsum(1); % Choose first element
                    [idx(1),idx(2),idx(3),idx(4),idx(5)] = ind2sub(size(modelPost), I_post(I_cumsum));
                    for iParam = 1:NumParams
                        simParams(jj, iParam) = params2fit(iParam, idx(iParam)); 
                    end
                    if strcmp(models{j}, 'fixed')
                        dataSim = changeprob_fixed_simulate(data, task(i), models{j}, simParams(jj,:));
                    elseif strcmp(models{j}, 'idealBayesian') || strcmp(models{j}, 'subBayesian_rlprior') || strcmp(models{j}, 'subBayesian_conservative') || strcmp(models{j}, 'subBayesian_pVec')
                        dataSim = changeprob_bocpd_simulate(data, task(i), models{j}, simParams(jj,:));
                    elseif strcmp(models{j}, 'exponential') || strcmp(models{j}, 'exponential_conservative')
                        dataSim = changeprob_exp_simulate(data, task(i), models{j}, simParams(jj,:));
                    elseif strcmp(models{j}, 'RL_probability') || strcmp(models{j}, 'RL_probability_conservative')
                        dataSim = changeprob_RLprob_simulate(data, task(i), models{j}, simParams(jj,:));
                    else
                        dataSim = changeprob_RLcriterion_simulate(data, task(i), models{j}, simParams(jj,:));
                    end
                    dataSim.SimParameters = simParams;
                    cd('/Users/elysenorton/Desktop/ChangeProb/matlab/ModelSimulations');
                    save(char(strcat('ChangeProb_Sim_', models{j}, '_', currentTask, '_', subID{ii}, '_', num2str(jj))), 'dataSim');
                    cd(currentFolder);
                end
            end
        end
    else
        for j = 1:numel(models)
            if j < 3
                parameters = [0 1 0 0 0 0 0 0];
            elseif j == 3 || j == 4 || j == 7
                parameters = [0 1 0 0 1 0 0 0];
            elseif j == 8
                parameters = [0 1 0 0 0 0 1 0];
            elseif j == 9
                parameters = [0 1 0 0 0 1 0 0];
            elseif j == 10
                parameters = [0 1 0 0 0 0 0 1];
            else
                parameters = [0 1 0 0 1 1 0 0];
            end
            NumParams = sum(parameters);
            I_params = find(parameters ~= 0);
            paramBounds = paramBounds_def(I_params,:);
            for iParam = 1:NumParams
                switch I_params(iParam)
                    case {1, 2}
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize); % sigma_noise
                    case 3
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize); % lambda
                    case 4
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize); % gamma
                    case 5
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize); % alpha
                    case 6
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize); % w
                    case 7
                        params2fit(iParam,:) = round(linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize)); % Tmax (discrete)
                    case 8
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize); % Range for minimum probability (pVec(1))
                end
            end
            for ii = 1:numel(subID_mixed)
                % Load data
                cd('/Users/elysenorton/Desktop/ChangeProb/data');
                load(strcat('ChangingProbabilitiesMixed_', subID_mixed{ii}));
                % Load fit data
                cd('/Users/elysenorton/Desktop/ChangeProb/matlab/ModelFit_data/SMALLlapse');
                load(char(strcat(subID_mixed{ii}, '_', models{j}, '_', currentTask)));
                cd(currentFolder);
                % Choose parameter vector with probability proportional to the joint posterior over parameters for given subject
                postVector = modelPost(:)/sum(modelPost(:));
                % Sort posterior vector
                [postSorted, I_post] = sortrows(postVector);
                % Compute the cumulative sum of sorted vector
                cumPost = cumsum(postSorted);
                for jj = 1:numSims
                    dataSim.randomSeed(jj,:) = i*j*ii*jj;
                    rng(dataSim.randomSeed(jj,:));
                    % Sample from posterior
                    I_cumsum = find(rand(1,1) < cumPost);
                    I_cumsum = I_cumsum(1); % Choose first element
                    [idx(1),idx(2),idx(3),idx(4),idx(5)] = ind2sub(size(modelPost), I_post(I_cumsum));
                    for iParam = 1:NumParams
                        simParams(jj, iParam) = params2fit(iParam, idx(iParam)); 
                    end
                    if strcmp(models{j}, 'fixed')
                        dataSim = changeprob_fixed_simulate(data, task(i), models{j}, simParams(jj,:));
                    elseif strcmp(models{j}, 'idealBayesian') || strcmp(models{j}, 'subBayesian_rlprior') || strcmp(models{j}, 'subBayesian_conservative') || strcmp(models{j}, 'subBayesian_pVec')
                        dataSim = changeprob_bocpd_simulate(data, task(i), models{j}, simParams(jj,:));
                    elseif strcmp(models{j}, 'exponential') || strcmp(models{j}, 'exponential_conservative')
                        dataSim = changeprob_exp_simulate(data, task(i), models{j}, simParams(jj,:));
                    elseif strcmp(models{j}, 'RL_probability') || strcmp(models{j}, 'RL_probability_conservative')
                        dataSim = changeprob_RLprob_simulate(data, task(i), models{j}, simParams(jj,:));
                    else
                        dataSim = changeprob_RLcriterion_simulate(data, task(i), models{j}, simParams(jj,:));
                    end
                    dataSim.SimParameters = simParams;
                    cd('/Users/elysenorton/Desktop/ChangeProb/matlab/ModelSimulations');
                    save(char(strcat('ChangeProb_Sim_', models{j}, '_', currentTask, '_', subID_mixed{ii}, '_', num2str(jj))), 'dataSim');
                    cd(currentFolder);
                end
            end
        end
    end
end

end

