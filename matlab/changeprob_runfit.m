function [marginalLikelihood, modelPost, nLL, rmse, fitParams, resp_model,...
    resp_obs, p_true, p_estimate, post] = changeprob_runfit(jobNumber)
%RUNFIT Runs model comparison for changing probability experiment
%   Detailed explanation goes here

    % Author:   Elyse norton
    % Email:    elyse.norton@gmail.com
    % Date:     11/5/2016

% Job number ranges from 1 to 154 (11 subjects x 2 tasks x 7 models)
if or(jobNumber < 1, jobNumber > (11*2*7))
    error('Please specify a number between 1 and 154');
end

subID = {'CWG', 'EGC', 'EHN', 'ERK', 'GK', 'HHL', 'JKT', 'JYZ', 'RND', 'SML', 'SQC'};
models = {'fixed', 'idealBayesian', 'exponential', 'RL_probability', ...
    'exponential_conservative', 'RL_probability_conservative', 'RL_criterion'};

% Which task?
if and(rem(jobNumber, 14) > 0, rem(jobNumber, 14) <= 7)
    task = 1;
    taskName = 'Overt';
else
    task = 2;
    taskName = 'Covert';
end

% Which model?
modelIndex = rem(jobNumber, 7);
if modelIndex == 0
    modelIndex = 7;
end
runModel = models{modelIndex};

% Which subject?
subIndex = ceil(jobNumber/14);
runSubject = subID{subIndex};

% Save file name
SaveFileName = strcat(runSubject, '_', runModel, '_', taskName);

% Fit data for subject, model, and task specified
parameters = [];

% Add project directory and subdirs to path
matlabdir = fileparts(which('changeprob_mL'));
basedir = matlabdir(1:find(matlabdir == filesep(), 1, 'last')-1);
addpath(genpath(basedir));
load(['ChangingProbabilities_', runSubject]); % Load data

if strcmp(runModel, 'exponential_conservative')
    runModel = 'exponential';
    if isempty(parameters)
        if task == 1
            parameters = [0 1 0 0 1 1];
        else
            parameters = [1 0 0 0 1 1];
        end
    end
elseif strcmp(runModel, 'RL_probability_conservative')
    runModel = 'RL_probability';
    if isempty(parameters)
        if task == 1
            parameters = [0 1 0 0 1 1];
        else
            parameters = [1 0 0 0 1 1];
        end
    end
end
[marginalLikelihood, modelPost, nLL, rmse, fitParams, resp_model,...
    resp_obs, p_true, p_estimate, post] = changeprob_mL(runModel, data, task, parameters);

save(SaveFileName, 'marginalLikelihood', 'modelPost', 'nLL', 'rmse', 'fitParams', ...
    'resp_model', 'resp_obs', 'p_true', 'p_estimate', 'post', ...
    'runSubject', 'runModel', 'subID', 'subIndex', 'taskName');

end

