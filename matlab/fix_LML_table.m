function [LML_table_correct] = fix_LML_table(LML_table, models)
%FIX_LML_TABLE Adds a correction to the original log marginal likelihood 
%scores for specified models and task
% Input:
    % LML_table - N x M table of log marginal likelihood scores where N is
    % the number of participants and M is the number of models
    % models - cell array of models to correct
% Output:
    % LML_table_correct - N x M table of corrected log marginal likelihood
    % scores
    
% Author:   Elyse Norton
% Email:    elyse.norton@gmail.com
% Date:     April/20/2019

models_default = {'idealBayesian', 'subBayesian_rlprior', 'subBayesian_pVec', ...
    'subBayesian_betahyp', 'fixed', 'exponential', 'exponential_conservative', ...
    'gold_nu', 'RL_criterion', 'behrens', 'behrens_conservative'};

% If no input variables are specified
if nargin < 1
    print('You must specify the LML_table to apply the correction');
end

% If no models are specified apply correction to all
if nargin < 2 || isempty(models)
    models = models_default;
end

% apply correction
LML_table_correct = LML_table - 2.1771;

% find the model indices
for i = 1:numel(models)
    idx_model(i) = find(strcmp(models{i}, models_default));
end

% select columns to return
LML_table_correct = LML_table_correct(:,idx_model);

end