function [lmL_Overt,lmL_Covert,initials,models,task] = changeprob_readFit(fitType,initials,models,task,dataPath)
%% CHANGEPROB_READFIT Read model fits for each observer, model, and task specified
    
    % Input
        % fitType: 'logmarglike' (default) or 'maxlike'
        % initials (Nx1): Observer identifier
        % models (Mx1): list of models to compare - relative model comparison scores
        % are computed relative to the first model in the list
        % task: overt (1) or covert (2) - compare performance in either or both tasks
        % dataPath: path of directory in which the data is stored - default
        % is current directory

    % Output
        % lmL_Overt (NxM): log marginal likelihood score for each observer and model
        % lmL_Covert (NxM): log marginal likelihood scores for each observer and model

%% Check inputs

    if nargin < 1 || isempty(fitType)   % Default - read LML files
        fitType = 'logmarglike';
    end
    
    switch fitType
        case 'logmarglike'; fileSuffix = [];
        case 'maxlike'; fileSuffix = '_maxlike';
        otherwise
            error('Unknown fitType - only ''logmarglike'' and ''maxlike'' are available.');
    end

    if nargin < 2 || isempty(initials) % Default - run all subjects
        initials = {'CWG', 'EGC', 'EHN', 'ERK', 'GK', 'HHL', 'JKT', 'JYZ', 'RND', 'SML', 'SQC'};
    end

    numInitials = numel(initials);

    if nargin < 3 || isempty(models) % Default - run all models
        models = {'fixed', 'idealBayesian', 'exponential', 'RL_probability', ...
            'exponential_conservative', 'RL_probability_conservative', 'RL_criterion', ...
            'subBayesian_rlprior', 'subBayesian_conservative', 'subBayesian_pVec', 'subBayesian_betahyp', ...
            'subBayesian_3param', 'gold', 'gold_nu', 'subBayesian_flex', 'behrens', 'behrens_conservative', ...
            'behrens_jump'};        
%         models = {'fixed', 'idealBayesian', 'exponential', 'RL_probability', ...
%             'exponential_conservative', 'RL_probability_conservative', 'RL_criterion', ...
%             'subBayesian_rlprior'};
    end

    numModels = numel(models);

    if nargin < 4 || isempty(task) % Default - run both tasks
        task = [1 2];
    end

    numTasks = numel(task);
    taskNames = {'Overt','Covert'};
    
    if nargin < 5 || isempty(dataPath)
        dataPath = pwd;
    end
    
    currentDir = pwd;

    % Plot model fits vs. observer data
    lmL_Overt = NaN(numInitials, numModels);
    lmL_Covert = lmL_Overt;
    
    for ii = 1:numTasks        
        for i = 1:numModels
            for j = 1:numInitials
                cd(dataPath); 
                filename = [initials{j}, '_', models{i}, '_', taskNames{task(ii)} fileSuffix '.mat'];
                if exist(filename,'file')
                    fit = load(filename);
                    switch task(ii)
                        case 1; lmL_Overt(j,i) = fit.logmargLikelihood;
                        case 2; lmL_Covert(j,i) = fit.logmargLikelihood;
                    end
                end
            end
        end        
    end

    cd(currentDir);
end
