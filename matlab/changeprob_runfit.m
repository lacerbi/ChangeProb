function [logmargLikelihood, modelPost, nLL, rmse, fitParams, resp_model,...
    resp_obs, p_true, p_estimate, post] = changeprob_runfit(jobNumber, fitType, fixNoise, gridSize, Nopts)
%RUNFIT Runs model comparison for changing probability experiment
%   Input:
    % jobNumber: unique model, task, data type, and subject identifier
    % fitType: fits data via log marginal
    % likelihood ('logmarglike') or max likelihood ('maxlike')
    % fixNoise: fixes the noise parameter
    % gridSize: determines the size of the grid for computing the log
    % marginal likelihood ('logmarglike' only)
    % Nopts: # optimization runs for maximum-likelihood fits ('maxlike' only)
    
% Output:
   % logmarglikelihood: log marginal likelihood
   % modelPost: the model posterior approximated over a grid and used to
   % compute the log marginal likelihood
   % nLL: negative log likelihood for the best fit parameters
   % rmse: root mean squared error between the model and data
   % fitParams: vector of the best fit parameters
   % resp_model: model prediction on each trial
   % resp_obs: observers response on each trial
   % p_true: true probability of category A on each trial
   % p_estimate: model estimate of the probability of category A on each
   % trial
   % post: ideal Bayesian model posterior

% Author:   Elyse norton
% Email:    elyse.norton@gmail.com
% Date:     10/3/2017

if nargin < 2; fitType = 'logmarglike'; end
if nargin < 3; fixNoise = []; end
if nargin < 4; gridSize = []; end
if nargin < 5; Nopts = []; end

subID = {'CWG', 'EGC', 'EHN', 'ERK', 'GK', 'HHL', 'JKT', 'JYZ', 'RND', 'SML', 'SQC'};
subID_mixed = {'CWG', 'EGC', 'EHN', 'ERK', 'HHL', 'RND', 'SML'}; % 7 of the 11 subjects also completed the mixed design experiment
models = {'fixed', 'idealBayesian', 'exponential', 'RL_probability', ...
    'exponential_conservative', 'RL_probability_conservative', 'RL_criterion', ...
    'subBayesian_rlprior', 'subBayesian_conservative', 'subBayesian_pVec', 'subBayesian_betahyp', ...
    'subBayesian_3param', 'gold', 'gold_nu', 'subBayesian_flex', 'behrens', 'behrens_conservative', ...
    'behrens_jump'}; % Models to fit
simModels = models; % Model used to simulate data

Nsubjs = numel(subID);
Nsubjs_mixed = numel(subID_mixed);
Nmodels = numel(models);
Ntasks = 2;     % Overt and covert
Nsims = 30;     % Number of simulations run for each model

% Job number ranges from 1 to 15,283 (11 subjects x 2 tasks x 17 models + 7 subjects x 17 models = 374 + 119 + 14,790)
    % Note: for each jobNumber > 493 we will fit all models to a complete
    % simulated data set at once
NdataJobs = (Nsubjs*Ntasks*Nmodels + Nsubjs_mixed*Nmodels);
NsimJobs = (Nsubjs*Ntasks*Nmodels + Nsubjs_mixed*Nmodels)*Nsims;
maxID = NdataJobs + NsimJobs;

if jobNumber < 1 || jobNumber > maxID
    error(['Please specify a number between 1 and ' num2str(maxID) '.']);
end

rng(jobNumber);     % Fix random seed

% Which subject?
if jobNumber <= Nsubjs*Ntasks*Nmodels 
    subIndex = rem(jobNumber-1,Nsubjs)+1;
    runSubject = subID{subIndex};
    simulatedData = [];
elseif jobNumber > Nsubjs*Ntasks*Nmodels && jobNumber <= NdataJobs
    jobNumber_mixed = jobNumber - Nsubjs*Ntasks*Nmodels;
    subIndex = rem(jobNumber_mixed-1,Nsubjs_mixed)+1;
    runSubject = subID_mixed{subIndex};
    simulatedData = [];
elseif jobNumber > NdataJobs
    simulatedData = 1;
    jobNumber_sim = jobNumber - NdataJobs;
    if jobNumber_sim <= Nsubjs*Ntasks*Nmodels*Nsims % Number between 1 and 9240
        % What is the sample number (1 to 30)?
        simNumIndex = ceil(jobNumber_sim/(Nsubjs*Ntasks*Nmodels));
        runSimNum = num2str(simNumIndex);
        % Update job number to be between 1 and 308
        jobNumber = mod(jobNumber_sim-1, Nsubjs*Ntasks*Nmodels)+1;
        % Which subject?
        subIndex = rem(jobNumber-1,Nsubjs)+1;
        runSubject = subID{subIndex};
    else
        jobNumber_mixed_sim = jobNumber_sim - Nsubjs*Ntasks*Nmodels*Nsims; % Number between 1 and 2940
        % What is the sample number (1 to 30)?
        simNumIndex = ceil(jobNumber_mixed_sim/(Nsubjs_mixed*Nmodels));
        runSimNum = num2str(simNumIndex);
        % Update job number to be between 1 and 70
        jobNumber_mixed = mod(jobNumber_mixed_sim-1, Nsubjs_mixed*Nmodels)+1;
        % Which subject?
        subIndex = rem(jobNumber_mixed-1,Nsubjs_mixed)+1;
        runSubject = subID_mixed{subIndex};
    end
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

% Simulation or fitting model?
if isempty(simulatedData)
    runModel = models{modelIndex};
else
    runSimModel = models{modelIndex};
end

% Add project directory and subdirs to path
matlabdir = fileparts(which('changeprob_logmarglike'));
basedir = matlabdir(1:find(matlabdir == filesep(), 1, 'last')-1);
addpath(genpath(basedir));
if isempty(simulatedData)
    if task == 3
        load(['ChangingProbabilitiesMixed_', runSubject]);
    else
        load(['ChangingProbabilities_', runSubject]); % Load data
    end
else
    load(['ChangeProb_Sim_', runSimModel, '_', taskName, '_', runSubject, '_', runSimNum]);
    data = dataSim;
    if isempty(gridSize)
        gridSize = 50;
    end
    simParams = data.SimParameters;
end

if isempty(simulatedData)
    initialRunModel = 1;
    NumRunModel = 1-initialRunModel;
else
%     initialRunModel = 1;
%     NumRunModel = Nmodels-initialRunModel;
    if strcmp(taskName, 'Mixed')
        if jobNumber_mixed < 71
            initialRunModel = 11;
            NumRunModel = Nmodels-initialRunModel;
        else
            initialRunModel = 1;
            NumRunModel = Nmodels-initialRunModel;
        end
    else
        if jobNumber < 221
            initialRunModel = 11;
            NumRunModel = Nmodels-initialRunModel;
        else
            initialRunModel = 1;
            NumRunModel = Nmodels-initialRunModel;
        end
    end
end

for ii = initialRunModel:(initialRunModel+NumRunModel)
    if isempty(simulatedData)
        if strcmp(fitType, 'maxlike')
            SaveFileName = strcat(runSubject, '_', runModel, '_', taskName, '_', fitType);
        else
            SaveFileName = strcat(runSubject, '_', runModel, '_', taskName);
        end
    else
        runModel = models{ii};
        SaveFileName = strcat(runSubject, '_', runSimModel, '_', runModel, '_', taskName, '_', runSimNum);
                
        try
            temp = load(SaveFileName);
            if ~isfield(temp,'simParams'); error('File does not contain simulated paramaters.'); end            
            % Fit already exists, skip to next fit
            fprintf('\n\n%s\n%s: Fit %d/%d already on file. Skipping.\n%s\n\n', repmat('#',[1,80]), datestr(now), ii, NumRunModel, repmat('#',[1,80]));
            continue;
        catch
            startTime = tic;
            fprintf('\n\n%s\n%s: Starting fit %d/%d.\n%s\n\n', repmat('#',[1,80]), datestr(now), ii, NumRunModel, repmat('#',[1,80]));
            % Either file does not exist or is corrupted
        end
        
    end
    parameters = [];
    switch(runModel)
        case 'exponential_conservative'
            runModel = 'exponential';
            if isempty(parameters)
                if task == 1 || task == 3
                    parameters = [0 1 0 0 1 1 0 0 0 0 0 0 0];
                else
                    parameters = [1 0 0 0 1 1 0 0 0 0 0 0 0];
                end
            end
        case 'RL_probability_conservative'
            runModel = 'RL_probability';
            if isempty(parameters)
                if task == 1 || task == 3
                    parameters = [0 1 0 0 1 1 0 0 0 0 0 0 0];
                else
                    parameters = [1 0 0 0 1 1 0 0 0 0 0 0 0];
                end
            end
        case 'subBayesian_rlprior'
            runModel = 'idealBayesian';
            if task == 1 || task == 3
                parameters = [0 1 0 0 0 0 1 0 0 0 0 0 0];
            else
                parameters = [1 0 0 0 0 0 1 0 0 0 0 0 0];
            end
        case 'subBayesian_conservative'
            runModel = 'idealBayesian';
            if task == 1 || task == 3
                parameters = [0 1 0 0 0 1 0 0 0 0 0 0 0];
            else
                parameters = [1 0 0 0 0 1 0 0 0 0 0 0 0];
            end
        case 'subBayesian_pVec'
            runModel = 'idealBayesian';
            if task == 1 || task == 3
                parameters = [0 1 0 0 0 0 0 1 0 0 0 0 0];
            else
                parameters = [1 0 0 0 0 0 0 1 0 0 0 0 0];
            end
        case 'subBayesian_betahyp'
            runModel = 'idealBayesian';
            if task == 1 || task == 3
                parameters = [0 1 0 0 0 0 0 0 1 0 0 0 0];
            else
                parameters = [1 0 0 0 0 0 0 0 1 0 0 0 0];
            end
        case 'subBayesian_3param'
            runModel = 'idealBayesian';
            if task == 1 || task == 3
                parameters = [0 1 0 0 0 0 1 0 1 0 0 0 0];
            else
                parameters = [1 0 0 0 0 0 1 0 1 0 0 0 0];
            end
            if isempty(gridSize)
                gridSize = 50;
            end
        case 'gold_nu'
            runModel = 'gold';
            if task == 1 || task == 3
                parameters = [0 1 0 0 0 0 0 0 0 1 1 0 1];
            else
                parameters = [1 0 0 0 0 0 0 0 0 1 1 0 1];
            end
            if isempty(gridSize)
                gridSize = 50;
            end
        case 'subBayesian_flex'
            runModel = 'idealBayesian';
            if task == 1 || task == 3
                parameters = [0 1 0 0 0 0 1 1 1 0 0 0 0 1 1];
            else
                parameters = [1 0 0 0 0 0 1 1 1 0 0 0 0 1 1];
            end
        case 'behrens_conservative'
            runModel = 'behrens_jump';
            if task == 1 || task == 3
                parameters = [0 1 0 0 0 1 0 0 0 0 0 0 0];
            else
                parameters = [1 0 0 0 0 1 0 0 0 0 0 0 0];
            end
            if isempty(gridSize)
                gridSize = 50;
            end
    end
    
    % Display model name, subject and task
    fprintf('Fitting model: %s. Subject: %s (#%d). Task: %s.\n\n',runModel,runSubject,subIndex,taskName);

    if strcmp(fitType, 'maxlike')
        [nLL, fitParams, resp_model, resp_obs, p_true, p_estimate,...
            logmargLikelihood,vbmc_fit] = changeprob_maxlike(runModel, data, task, parameters, [], Nopts);
        rmse = [];
        post = [];
        modelPost = [];
        v_estimate = [];
    else
        [logmargLikelihood, modelPost, nLL, rmse, fitParams, resp_model,...
            resp_obs, p_true, p_estimate, post, v_estimate] = changeprob_logmarglike(runModel, data, task, parameters, gridSize, [], fixNoise, simulatedData);
        vbmc_fit = [];
    end

    fprintf('MAP parameters:\n');
    fitParams

    if isempty(simulatedData)
        save(SaveFileName, 'logmargLikelihood', 'modelPost', 'nLL', 'rmse', 'fitParams', ...
            'resp_model', 'resp_obs', 'p_true', 'p_estimate', 'post', 'v_estimate', ...
            'runSubject', 'runModel', 'subID', 'subIndex', 'taskName','vbmc_fit');
    else
        save(SaveFileName, 'logmargLikelihood', 'modelPost', 'nLL', 'rmse', 'fitParams', ...
            'resp_model', 'resp_obs', 'p_true', 'p_estimate', 'post', 'v_estimate', ...
            'runSubject', 'runModel', 'subID', 'subIndex', 'taskName', 'vbmc_fit', 'runSimNum', 'runSimModel', 'simParams');
        fprintf('\n\n%s\n%s: Finished fit %d/%d (running time = %.2f h).\n%s\n\n', repmat('#',[1,80]), datestr(now), ii, NumRunModel, toc(startTime)/3600, repmat('#',[1,80]));       
    end
end

end

