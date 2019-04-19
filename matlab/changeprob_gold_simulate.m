function dataSim = changeprob_gold_simulate(data, task, parameters)
%CHANGEPROB_GOLD_SIMULATE Simulates an observer with a mixture of delta rules (Wilson et al., 2013).
% (Documentation to be written.)
%
% Author:   Elyse Norton
% Email:    elyse.norton@gmail.com
% Date:     August 2017

% Check input arguments

% Check for data file
if nargin < 1; data = []; task = 1; parameters = []; end

% Determine task to simulate
if nargin < 2; task = 1; parameters = []; end
if isempty(task); task = 1; end
if task ~= 1 && task ~= 2 && task ~=3; error('TASK can only be 1 (overt-criterion), 2 (covert-criterion), or 3 (mixed design).'); end

% Parameter vector: #1 is SIGMA_ELLIPSE / SIGMA_CRITERION, #2 and 3 are delta1 and delta2, and
% #4 is nu_p and #5 is hRate (optional)
if nargin < 3 || isempty(parameters); parameters = [5 1.5 4.5 log(2) .01]; end

% Get session parameters
if isempty(data) 
    [NumTrials, sigma_ellipseData, mu, sigma_s, C, S, p_true, ~, ~] = changeprob_getSessionParameters([], task, [parameters(1), parameters(1)]); % Generate fake data using specified parameters
else
    [NumTrials, sigma_ellipseData, mu, sigma_s, C, S, p_true, ~, ~] = changeprob_getSessionParameters(data, task); % Determine session parameters from provided data set
end

switch numel(parameters)
    case {1,2}
        error('Not enough input parameters. Model must have at least 3 parameters');
    case 3
        if task ~= 2
            sigma_criterion = parameters(1);
            sigma_ellipse = sigma_ellipseData;
        else
            sigma_ellipse = parameters(1);
        end
        delta1 = parameters(2)^2;
        delta2 = parameters(3)^2;
        nu_p = 2;
        hRate = .01; % Default hazard rate
    case 4
        if task ~= 2
            sigma_criterion = parameters(1);
            sigma_ellipse = sigma_ellipseData;
        else
            sigma_ellipse = parameters(1);
        end
        delta1 = parameters(2)^2;
        delta2 = parameters(3)^2;
        nu_p = exp(parameters(4));
        hRate = .01;
    case 5
        if task ~= 2
            sigma_criterion = parameters(1);
            sigma_ellipse = sigma_ellipseData;
        else
            sigma_ellipse = parameters(1);
        end
        delta1 = parameters(2)^2;
        delta2 = parameters(3)^2;
        nu_p = exp(parameters(4));
        hRate = parameters(5);
    otherwise
        error('Too many input parameters');
end
sigma = sqrt(sigma_s^2 + sigma_ellipse^2);
nodes = [1, 1+delta1, delta1+delta2+1];
numNodes = numel(nodes);

%% Start loop over trials

p_initial = .5;
z_model = zeros(1,NumTrials);
resp_obs = z_model;
score = z_model;
beta = z_model;
p_estNode = zeros(NumTrials, numNodes); % Probability estimate at each node (updated on each trial)
likelihood_dataNode = p_estNode; % Likelihood of the data given each node
changept_prior = zeros(numNodes, numNodes); % Change-point prior
p_nodeWeight = p_estNode; % Weights given to each node (updated on each trial)
p_estimate = zeros(NumTrials, 1); % Overall estimate of pA computed by taking a weighted average between the probability esitmates at each node

p_estNode(1,:) = p_estNode(1,:) + p_initial;
p_nodeWeight(1,:) = [1 0 0]; % Assume a change just occurred
p_estimate(1,:) = p_estNode(1,:)*p_nodeWeight(1,:)';
beta(1) = p_estimate(1,:)/(1-p_estimate(1,:));
z_model(:,1) = sigma^2 * log(beta(1)) / diff(mu);

if task == 1
    resp_obs(:,1) = z_model(:,1) + randn(1,1)*sigma_criterion;
    if and(C(1) == 1, S(1) <= resp_obs(:,1)) || and(C(1) == 0, S(1) > resp_obs(:,1))
        score(:,1) = 1;
    else
        score(:,1) = 0;
    end
else
    if S(1)+randn(1,1)*sigma_ellipse <= z_model(:,1)
        resp_obs(:,1) = 1;
    else
        resp_obs(:,1) = 0;
    end
    score(:,1) = double(resp_obs(:,1) == C(1));
end

Cprev = C(1);

for t = 2:NumTrials
    for idx_node = 1:numNodes % Subscript i in text - end node state
        % Update the probability estimate at each node
        p_estNode(t,idx_node) = p_estNode(t-1,idx_node) + (1/(nodes(idx_node)+nu_p))*(Cprev-p_estNode(t-1,idx_node));
        % Compute the likelihood of the new data given each node
        likelihood_dataNode(t,idx_node) = p_estNode(t,idx_node)^(Cprev==1)*(1-p_estNode(t,idx_node))^(Cprev~=1);
        % Compute the change-point prior
        for idx_startNode = 1:numNodes % Subscript j in text - start node state
            % Change-point prior: Change
            if idx_node == 1
                p_change = 1; % Probability of a change (only equals one if idx_node = 1, otherwise equals 0)
            else
                p_change = 0;
            end
            % Change-point prior: No change (advance to next node or stay put)
            if idx_startNode ~= numNodes
                if idx_node == idx_startNode
                    p_noChange = (nodes(idx_startNode+1) - nodes(idx_startNode) - 1)/(nodes(idx_startNode+1)-nodes(idx_startNode));
                elseif idx_node == idx_startNode + 1
                    p_noChange = 1/(nodes(idx_startNode+1)-nodes(idx_startNode));
                else
                    p_noChange = 0;
                end
            elseif idx_startNode == numNodes
                if idx_node == numNodes
                    p_noChange = 1; % Self-transition probability at final node (can't advance)
                else
                    p_noChange = 0;
                end
            end
            % Take a weighted average of the change and no change
            % probabilities, weighted by h and (1-h), respectively
            changept_prior(idx_node, idx_startNode) = hRate*p_change + (1-hRate)*p_noChange;
        end
    end
    % Update weight for each node's probability estimate
    p_nodeWeight(t,:) = likelihood_dataNode(t,:)'.*(changept_prior*p_nodeWeight(t-1,:)');
    % Normalize weights so they sum to 1
    p_nodeWeight(t,:) = p_nodeWeight(t,:)./sum(p_nodeWeight(t,:));
    % Take a weighted average of the probability estimates for each node
    p_estimate(t,:) = p_estNode(t,:)*p_nodeWeight(t,:)';
    % Used probability estimate to determine the criterion
    beta(t) = p_estimate(t,:)/(1-p_estimate(t,:));
    z_model(:,t) = sigma^2 * log(beta(t)) / diff(mu);
    % Response and score
    if task == 1 || and(task == 3, mod(t,5)==0)
        resp_obs(:,t) = z_model(:,t) + randn(1,1)*sigma_criterion;
        if and(C(t) == 1, S(t) <= resp_obs(:,t)) || and(C(t) == 0, S(t) > resp_obs(:,t))
            score(:,t) = 1;
        else
            score(:,t) = 0;
        end
    else
        if S(t)+randn(1,1)*sigma_ellipse <= z_model(:,t)
            resp_obs(:,t) = 1;
        end
        score(:,t) = double(resp_obs(:,t) == C(t));
    end
    Cprev = C(t);
end

dataSim.NumTrials = NumTrials;
dataSim.mu = mu;
dataSim.sigma_s = sigma_s;
dataSim.sigma = sigma;
dataSim.sigmaEllipse = sigma_ellipse;
dataSim.Category = C;
dataSim.Stimulus = S;
dataSim.pA = p_true;
dataSim.response = resp_obs;
dataSim.score = score;

end