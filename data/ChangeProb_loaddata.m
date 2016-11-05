function data = ChangeProb_loaddata()
%CHANGEPROB_LOADDATA Load all datasets from changeprobability experiment.

subID = {'CWG', 'EGC', 'EHN', 'ERK', 'GK', 'HHL', 'JKT', 'JYZ', 'RND', 'SML', 'SQC'};

data = [];
for i = 1:numel(subID)
    temp = load(['ChangingProbabilities_' subID{i} '.mat']);
    data{i} = temp.data;    
end

end