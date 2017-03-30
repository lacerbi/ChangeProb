function [nLL, rmse, p_estimate, resp_model, p_nodeWeight] = changeprob_gold_nll(parameters, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task)
%CHANGEPROB_GOLD_NLL "A Mixture of Delta-Rules Approximation to Bayesian
%Inference in Change-Point Problems" - Wilson et al. (2013). PLoS Computational Biology
% (Documentation to be written.)
%
% Author:   Elyse Norton
% Email:    elyse.norton@gmail.com
% Date:     Mar/27/2017

% Generate fake data and set default parameter vector
% Parameter vector: #1 is SIGMA_ELLIPSE, #2 is SIGMA_CRITERION, #3 is LAPSE, 
% #4 is GAMMA, #5 is ALPHA, and #6 is W, #7 is rl_prior, #8 is pVec, #9 is betahyp,
% #10 is delta node 1, #11 is delta node 2, #12 is hazard rate
if nargin < 1 || isempty(parameters)
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs, score] = changeprob_getSessionParameters();
    parameters = [sigma_ellipse, sigma_ellipse, 1e-4, Inf, .2, 1, 0, 0, 0, 2, 5, .01];
    task = 1; 
end

% Generate fake data using parameters specified
if nargin < 2
    task = 1;
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs, score] = changeprob_getSessionParameters([], task, [parameters(1), parameters(2)]);
end

% Make sure all the experiment variables are specified
if nargin > 1 && nargin < 9
    error('Missing input arguments. You must specify all or no experimental variables.');
end

% Task (1 overt, 2 covert)
if nargin == 9 || isempty(task); task = 1; end
if task ~= 1 && task ~= 2 && task ~=3; error('TASK can only be 1 (overt-criterion), 2 (covert-criterion), or 3 (mixed design).'); end

% Observer model parameters

sigma_ellipse = parameters(1);
sigma_criterion = parameters(2);
lambda = parameters(3);
delta1 = parameters(10)^2;
delta2 = parameters(11)^2;
nodes = [1, 1+delta1, delta1+delta2+1];
hRate = parameters(12);
nu_p = 2;
numNodes = numel(nodes);

%% Start loop over trials

p_initial = .5;
z_model = zeros(1,NumTrials);
PChatA = z_model;
log_P = z_model;
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
    log_P(:,1) = -0.5*log(2*pi*sigma_criterion) - 0.5*((resp_obs(1)-z_model(:,1))./sigma_criterion).^2;
else
    PChatA(:,1) = 1 - 0.5*erfc( -(S(1) - z_model(1)) / (sqrt(2)*sigma_ellipse) );
    log_P(:,1) = log(PChatA(:,1)).*(resp_obs(1)==1) + log(1-PChatA(:,1)).*(resp_obs(1)~=1);
end
Cprev = C(1);

for t = 2:NumTrials
    for idx_node = 1:numNodes
        % Update the probability estimate at each node
        p_estNode(t,idx_node) = p_estNode(t-1,idx_node) + (1/(nodes(idx_node)+nu_p))*(Cprev-p_estNode(t-1,idx_node));
        % Compute the likelihood of the new data given each node
        likelihood_dataNode(t,idx_node) = p_estNode(t,idx_node)^(Cprev==1)*(1-p_estNode(t,idx_node))^(Cprev~=1);
        % Compute the change-point prior
        endNodeState = nodes(idx_node); % Subscript i in text
        for idx_startNode = 1:numNodes % Subscript j in text
            startNodeState = nodes(idx_startNode);
            % Change-point prior: Change
            if idx_node == 1
                p_change = 1; % Probability of a change (only equals one if idx_node = 1, otherwise equals 0)
            else
                p_change = 0;
            end
            % Change-point prior: No change
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
                    p_noChange = 1; % Self-transition probability at final node
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
    % Compute the log likelihood
    if task == 1 || and(task == 3, mod(t,5)==0)
        log_P(:,t) = -0.5*log(2*pi*sigma_criterion) - 0.5*((resp_obs(t)-z_model(:,t))./sigma_criterion).^2;
        if lambda > 0
            log_P(:,t) = log(lambda/360 + (1-lambda)*exp(log_P(:,t)));
        end
    else
        PChatA(:,t) = 1 - 0.5*erfc( -(S(t) - z_model(t)) / (sqrt(2)*sigma_ellipse) );
        if lambda > 0
            PChatA(:,t) = lambda/2 + (1-lambda)*PChatA(:,t);
        end
        log_P(:,t) = log(PChatA(:,t)).*(resp_obs(t)==1) + log(1-PChatA(:,t)).*(resp_obs(t)~=1);
    end
    Cprev = C(t);
end

if task == 1
    resp_model = z_model;
elseif task ==2
    resp_model = PChatA;
else
    resp_model = PChatA;
    resp_model(5:5:NumTrials) = z_model(5:5:NumTrials);
end

nLL = -nansum(log_P);
rmse = sqrt(mean((p_estimate - p_true).^2));

end