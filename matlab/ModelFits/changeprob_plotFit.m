function [lmL_Overt, lmL_Covert] = changeprob_plotFit(initials, models, task, dataPath, save)
%% CHANGEPROB_PLOTFIT Plots the model fits for each observer, model, and task specified
    % When more than one model is specified it also plots the relative
    % model comparison scores (relative to the first model specified) in a
    % bar plot
    
    % Input
        % initials (Nx1): Observer identifier
        % models (Mx1): list of models to compare - relative model comparison scores
        % are computed relative to the first model in the list
        % task: overt (1) or covert (2) - compare performance in either or both tasks
        % dataPath: path of directory in which the data is stored - default
        % is current directory
        % save: path of directory you want to save the figures in

    % Output
        % lmL_Overt (NxM): log marginal likelihood score for each observer and model
        % lmL_Covert (NxM): log marginal likelihood scores for each observer and model

%% Check inputs

    if nargin < 1 || isempty(initials) % Default - run all subjects
        initials = {'CWG', 'EGC', 'EHN', 'ERK', 'GK', 'HHL', 'JKT', 'JYZ', 'RND', 'SML', 'SQC'};
    end

    numInitials = numel(initials);

    if nargin < 2 || isempty(models) % Default - run all models
        models = {'fixed', 'idealBayesian', 'exponential', 'RL_probability', ...
            'exponential_conservative', 'RL_probability_conservative', 'RL_criterion', ...
            'subBayesian_rlprior'};
    end

    numModels = numel(models);

    if nargin < 3 || isempty(task) % Default - run both tasks
        task = [1 2];
    end

    numTasks = numel(task);
    
    if nargin < 4 || isempty(dataPath)
        dataPath = pwd;
    end
    
    currentDir = pwd;
    
    if nargin < 5 || isempty(save) % Default - don't save
        save = [];
    end

    % Plot model fits vs. observer data
    lmL_Overt = zeros(numInitials, numModels);
    lmL_Covert = lmL_Overt;
    [R,C] = getRowsCols(numInitials,'R');
    h = 1;
    for ii = 1:numTasks
        if task(ii) == 1
            for i = 1:numModels
                f = figure(h);
                figName = char(strcat(models{i}, {' '}, 'Overt'));
                set(f,'name',figName,'numbertitle','off')
                for j = 1:numInitials
                    cd(dataPath);
                    load(strcat(initials{j}, '_', models{i}, '_', 'Overt'));
                    lmL_Overt(j,i) = logmargLikelihood;
                    subplot(R,C,j);
                    hold on; title(initials{j}); 
                    xlabel('Trial number'); 
                    ylabel('Orientation (deg)');
                    plot(1:800, smooth(resp_obs, 11), '-k', 'linewidth', 1);
                    plot(1:800, smooth(resp_model, 11), '-r', 'linewidth', 1);
                    hold off;
                end
                if ~isempty(save)
                    cd(save);
                    print(strcat(models{i}, '_', 'Overt'), '-dpdf');
                end
                h = h+1;
            end
        else
            for i = 1:numModels
                f = figure(h);
                set(f,'name',char(strcat(models{i}, {' '}, 'Covert')),'numbertitle','off');
                for j = 1:numInitials
                    cd(dataPath);
                    load(strcat(initials{j}, '_', models{i}, '_', 'Covert'));
                    lmL_Covert(j,i) = logmargLikelihood;
                    subplot(R,C,j);
                    hold on; title(initials{j}); 
                    xlabel('Trial number'); 
                    ylabel('P(A)');
                    axis([0 800 0 1]);
                    plot(1:800, smooth(resp_obs, 11), '-k', 'linewidth', 1);
                    plot(1:800, smooth(resp_model, 11), '-r', 'linewidth', 1);
                    hold off;
                end
                if ~isempty(save)
                    cd(save);
                    print(strcat(models{i}, '_', 'Covert'), '-dpdf');
                end
                h = h+1;
            end
        end
    end

    % Plot bar graph for model comparison if multiple models are called
    if numModels > 1 
        if numTasks == 2
            modelRelScore_Overt = lmL_Overt-lmL_Overt(:,1);
            modelRelScore_Covert = lmL_Covert-lmL_Covert(:,1);
            f = figure(h); 
            set(f,'name','Relative model fits','numbertitle','off');
            subplot(2,1,1); hold on;
            xlabel('Observer');
            ylabel('Relative log marginal likelihood (Covert)');
            bar([modelRelScore_Covert(:,2:end); mean(modelRelScore_Covert(:,2:end))]);
            axis([0 numInitials+2 -10 150]);
            set(gca,'xTick',1:(numInitials+1));
            initials{end+1} = 'Avg';
            set(gca, 'xTickLabel', initials);
            subplot(2,1,2); hold on;
            xlabel('Observer');
            ylabel('Relative log marginal likelihood (Overt)');
            bar([modelRelScore_Overt(:,2:end); mean(modelRelScore_Overt(:,2:end))]);
            axis([0 13 -10 150]);
            set(gca,'xTick',1:12);
            set(gca, 'xTickLabel', initials);
            legend(models(2:end));
            if ~isempty(save)
                cd(save);
                print('LogMarginalLikelihood', '-dpdf');
            end
        elseif task == 2
            modelRelScore_Covert = lmL_Covert-lmL_Covert(:,1);
            f = figure(h); hold on;
            set(f,'name','Relative model fits','numbertitle','off');
            xlabel('Observer');
            ylabel('Relative log marginal likelihood (Covert)');
            bar([modelRelScore_Covert(:,2:end); mean(modelRelScore_Covert(:,2:end))]);
            axis([0 numInitials+2 -10 150]);
            set(gca,'xTick',1:(numInitials+1));
            initials{end+1} = 'Avg';
            set(gca, 'xTickLabel', initials);
            legend(models(2:end), 'location', 'best');
            if ~isempty(save)
                cd(save);
                print('LogMarginalLikelihood_Covert', '-dpdf');
            end
        elseif task == 1
            modelRelScore_Overt = lmL_Overt-lmL_Overt(:,1);
            f = figure(h); hold on;
            set(f,'name','Relative model fits','numbertitle','off');
            xlabel('Observer');
            ylabel('Relative log marginal likelihood (Overt)');
            bar([modelRelScore_Overt(:,2:end); mean(modelRelScore_Overt(:,2:end))]);
            axis([0 numInitials+2 -10 150]);
            set(gca,'xTick',1:(numInitials+1));
            initials{end+1} = 'Avg';
            set(gca, 'xTickLabel', initials);
            legend(models(2:end));
            if ~isempty(save)
                cd(save);
                print('LogMarginalLikelihood_Overt', '-dpdf');
            end
        end
    end
    cd(currentDir);
end




