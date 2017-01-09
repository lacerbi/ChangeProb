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
            for ii = 1:numel(subID)
                % Load data
                cd('/Users/elysenorton/Desktop/ChangeProb/data');
                load(strcat('ChangingProbabilities_', subID{ii}));
                % Load fit data
                cd('/Users/elysenorton/Desktop/ChangeProb/matlab/ModelFit_data/SMALLlapse');
                load(char(strcat(subID{ii}, '_', models{j}, '_', currentTask)));
                cd(currentFolder);
                for jj = 1:numSims
                    if strcmp(models{j}, 'fixed')
                        dataSim = changeprob_fixed_simulate(data, task(i), models{j}, fitParams);
                    elseif strcmp(models{j}, 'idealBayesian') || strcmp(models{j}, 'subBayesian_rlprior') || strcmp(models{j}, 'subBayesian_conservative') || strcmp(models{j}, 'subBayesian_pVec')
                        dataSim = changeprob_bocpd_simulate(data, task(i), models{j}, fitParams);
                    elseif strcmp(models{j}, 'exponential') || strcmp(models{j}, 'exponential_conservative')
                        dataSim = changeprob_exp_simulate(data, task(i), models{j}, fitParams);
                    elseif strcmp(models{j}, 'RL_probability') || strcmp(models{j}, 'RL_probability_conservative')
                        dataSim = changeprob_RLprob_simulate(data, task(i), models{j}, fitParams);
                    else
                        dataSim = changeprob_RLcriterion_simulate(data, task(i), models{j}, fitParams);
                    end
                    cd('/Users/elysenorton/Desktop/ChangeProb/matlab/ModelSimulations');
                    save(char(strcat('ChangeProb_Sim_', models{j}, '_', currentTask, '_', subID{ii}, '_', num2str(jj))), 'dataSim');
                    cd(currentFolder);
                end
            end
        end
    else
        for j = 1:numel(models)
            for ii = 1:numel(subID_mixed)
                % Load data
                cd('/Users/elysenorton/Desktop/ChangeProb/data');
                load(strcat('ChangingProbabilitiesMixed_', subID_mixed{ii}));
                % Load fit data
                cd('/Users/elysenorton/Desktop/ChangeProb/matlab/ModelFit_data/SMALLlapse');
                load(char(strcat(subID_mixed{ii}, '_', models{j}, '_', currentTask)));
                cd(currentFolder);
                for jj = 1:numSims
                    if strcmp(models{j}, 'fixed')
                        dataSim = changeprob_fixed_simulate(data, task(i), models{j}, fitParams);
                    elseif strcmp(models{j}, 'idealBayesian') || strcmp(models{j}, 'subBayesian_rlprior') || strcmp(models{j}, 'subBayesian_conservative') || strcmp(models{j}, 'subBayesian_pVec')
                        dataSim = changeprob_bocpd_simulate(data, task(i), models{j}, fitParams);
                    elseif strcmp(models{j}, 'exponential') || strcmp(models{j}, 'exponential_conservative')
                        dataSim = changeprob_exp_simulate(data, task(i), models{j}, fitParams);
                    elseif strcmp(models{j}, 'RL_probability') || strcmp(models{j}, 'RL_probability_conservative')
                        dataSim = changeprob_RLprob_simulate(data, task(i), models{j}, fitParams);
                    else
                        dataSim = changeprob_RLcriterion_simulate(data, task(i), models{j}, fitParams);
                    end
                    cd('/Users/elysenorton/Desktop/ChangeProb/matlab/ModelSimulations');
                    save(char(strcat('ChangeProb_Sim_', models{j}, '_', currentTask, '_', subID_mixed{ii}, '_', num2str(jj))), 'dataSim');
                    cd(currentFolder);
                end
            end
        end
    end
end

end

