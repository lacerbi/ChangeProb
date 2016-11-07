function [marginalLikelihood, modelPost, nLL, rmse, fitParams, resp_model,...
    resp_obs, p_true, p_estimate, post] = changeprob_runfit(jobNumber)
%RUNFIT Runs model comparison for changing probability experiment
%   Detailed explanation goes here

% Author:   Elyse norton
% Email:    elyse.norton@gmail.com
% Date:     11/5/2016

% Submit RL models
% ./submitfit.sh 2 5,19,33,47,61,75,89,103,117,131,145
% ./submitfit.sh 3 6,20,34,48,62,76,90,104,118,132,146

subID = {'CWG', 'EGC', 'EHN', 'ERK', 'GK', 'HHL', 'JKT', 'JYZ', 'RND', 'SML', 'SQC'};
models = {'fixed', 'idealBayesian', 'exponential', 'RL_probability', ...
    'exponential_conservative', 'RL_probability_conservative', 'RL_criterion'};

Nsubjs = numel(subID);
Nmodels = numel(models);
Ntasks = 2;     % Overt and covert

% Job number ranges from 1 to 154 (11 subjects x 2 tasks x 7 models)
maxID = Nsubjs*Ntasks*Nmodels;
if jobNumber < 1 || jobNumber > maxID
    error(['Please specify a number between 1 and ' num2str(maxID) '.']);
end

% Which subject?
subIndex = rem(jobNumber-1,Nsubjs)+1;
runSubject = subID{subIndex};

% Which task?
if rem(jobNumber-1, Nsubjs*Ntasks) < Nsubjs
    task = 1;
    taskName = 'Overt';
else
    task = 2;
    taskName = 'Covert';
end

% Which model?
modelIndex = ceil(jobNumber/(Nsubjs*Ntasks));
runModel = models{modelIndex};

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

[logmargLikelihood, modelPost, nLL, rmse, fitParams, resp_model,...
    resp_obs, p_true, p_estimate, post] = changeprob_logmarglike(runModel, data, task, parameters);

save(SaveFileName, 'logmargLikelihood', 'modelPost', 'nLL', 'rmse', 'fitParams', ...
    'resp_model', 'resp_obs', 'p_true', 'p_estimate', 'post', ...
    'runSubject', 'runModel', 'subID', 'subIndex', 'taskName');

end

