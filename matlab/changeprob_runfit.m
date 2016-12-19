function [logmargLikelihood, modelPost, nLL, rmse, fitParams, resp_model,...
    resp_obs, p_true, p_estimate, post] = changeprob_runfit(jobNumber, fixNoise, gridSize)
%RUNFIT Runs model comparison for changing probability experiment
%   Detailed explanation goes here

% Author:   Elyse norton
% Email:    elyse.norton@gmail.com
% Date:     12/15/2016

% Submit RL models
% ./submitfit.sh 2 5,19,33,47,61,75,89,103,117,131,145
% ./submitfit.sh 3 6,20,34,48,62,76,90,104,118,132,146

if nargin < 2; fixNoise = []; end
if nargin < 3; gridSize = []; end

subID = {'CWG', 'EGC', 'EHN', 'ERK', 'GK', 'HHL', 'JKT', 'JYZ', 'RND', 'SML', 'SQC'};
subID_mixed = {'CWG', 'EGC', 'EHN', 'ERK', 'HHL', 'RND', 'SML'}; % 7 of the 11 subjects also completed the mixed design experiment
models = {'fixed', 'idealBayesian', 'exponential', 'RL_probability', ...
    'exponential_conservative', 'RL_probability_conservative', 'RL_criterion', ...
    'subBayesian_rlprior', 'subBayesian_conservative', 'subBayesian_pVec'};

Nsubjs = numel(subID);
Nsubjs_mixed = numel(subID_mixed);
Nmodels = numel(models);
Ntasks = 2;     % Overt and covert

% Job number ranges from 1 to 290 (11 subjects x 2 tasks x 10 models + 7 subjects x 10 models)
maxID = Nsubjs*Ntasks*Nmodels + Nsubjs_mixed*Nmodels;
maxID_mixed = Nsubjs_mixed*Nmodels;
if jobNumber < 1 || jobNumber > maxID
    error(['Please specify a number between 1 and ' num2str(maxID) '.']);
end

rng(jobNumber);     % Fix random seed

% Which subject?
if jobNumber <= Nsubjs*Ntasks*Nmodels 
    subIndex = rem(jobNumber-1,Nsubjs)+1;
    runSubject = subID{subIndex};
else
    jobNumber_mixed = jobNumber - Nsubjs*Ntasks*Nmodels;
    subIndex = rem(jobNumber_mixed-1,Nsubjs_mixed)+1;
    runSubject = subID_mixed{subIndex};
end

% Which task?
if jobNumber <= Nsubjs*Ntasks*Nmodels
    if rem(jobNumber-1, Nsubjs*Ntasks) < Nsubjs
        task = 1;
        taskName = 'Overt';
    else
        task = 2;
        taskName = 'Covert';
    end
else
    task = 3;
    taskName = 'Mixed';
end

% Which model?
if jobNumber <= Nsubjs*Ntasks*Nmodels
    modelIndex = ceil(jobNumber/(Nsubjs*Ntasks));
else
    modelIndex = ceil(jobNumber_mixed/Nsubjs_mixed);
end
runModel = models{modelIndex};

% Save file name
SaveFileName = strcat(runSubject, '_', runModel, '_', taskName);

% Fit data for subject, model, and task specified
parameters = [];

% Add project directory and subdirs to path
matlabdir = fileparts(which('changeprob_mL'));
basedir = matlabdir(1:find(matlabdir == filesep(), 1, 'last')-1);
addpath(genpath(basedir));
if task == 3
    load(['ChangingProbabilitiesMixed_', runSubject]);
else
    load(['ChangingProbabilities_', runSubject]); % Load data
end

switch(runModel)
    case 'exponential_conservative'
        runModel = 'exponential';
        if isempty(parameters)
            if task == 1 || task == 3
                parameters = [0 1 0 0 1 1 0 0];
            else
                parameters = [1 0 0 0 1 1 0 0];
            end
        end
    case 'RL_probability_conservative'
        runModel = 'RL_probability';
        if isempty(parameters)
            if task == 1 || task == 3
                parameters = [0 1 0 0 1 1 0 0];
            else
                parameters = [1 0 0 0 1 1 0 0];
            end
        end
    case 'subBayesian_rlprior'
        runModel = 'idealBayesian';
        if task == 1 || task == 3
            parameters = [0 1 0 0 0 0 1 0];
        else
            parameters = [1 0 0 0 0 0 1 0];
        end
    case 'subBayesian_conservative'
        runModel = 'idealBayesian';
        if task == 1 || task == 3
            parameters = [0 1 0 0 0 1 0 0];
        else
            parameters = [1 0 0 0 0 1 0 0];
        end
    case 'subBayesian_pVec'
        runModel = 'idealBayesian';
        if task == 1 || task == 3
            parameters = [0 1 0 0 0 0 0 1];
        else
            parameters = [1 0 0 0 0 0 0 1];
        end
end

[logmargLikelihood, modelPost, nLL, rmse, fitParams, resp_model,...
    resp_obs, p_true, p_estimate, post] = changeprob_logmarglike(runModel, data, task, parameters, gridSize, [], fixNoise);

fprintf('MAP parameters:\n');
fitParams

save(SaveFileName, 'logmargLikelihood', 'modelPost', 'nLL', 'rmse', 'fitParams', ...
    'resp_model', 'resp_obs', 'p_true', 'p_estimate', 'post', ...
    'runSubject', 'runModel', 'subID', 'subIndex', 'taskName');

end

