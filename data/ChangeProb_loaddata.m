function data = ChangeProb_loaddata()
%CHANGEPROB_LOADDATA Load all datasets from changeprobability experiment.

subjects = {'EHN','ERK','GK','JKT','SML','EGC','HHL','JYZ','RND','SQC','CWG'};

data = [];
for i = 1:numel(subjects)
    temp = load(['ChangingProbabilities_' subjects{i} '.mat']);
    data{i} = temp.data;    
end

end