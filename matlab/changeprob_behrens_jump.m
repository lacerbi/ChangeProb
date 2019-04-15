function [nLL, rmse, p_estimate, resp_model, v_estimate] = changeprob_behrens_jump(parameters, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task)
%CHANGEPROB_BEHRENS_JUMP Bayesian observer model adapted from Behrens et al. (2007)
%paper titled "Learning the value of information in an uncertain world"
%assuming knowledge of jumps in probability
% (Documentation to be written.)
%
% Author:   Elyse Norton
% Email:    elyse.norton@gmail.com
% Date:     March/03/2019

% Generate fake data and set default parameter vector
% Parameter vector: #1 is SIGMA_ELLIPSE, #2 is SIGMA_CRITERION, #3 is LAPSE, 
% #4 is GAMMA (not used here), #5 is alpha (not used here), and #6 is W
if nargin < 1 || isempty(parameters)
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs] = changeprob_getSessionParameters();
    parameters = [sigma_ellipse, sigma_ellipse, 1e-4, Inf, .2, 1, 0];
    task = 1; 
end

% Generate fake data using parameters specified
if nargin < 2
    task = 1;
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs] = changeprob_getSessionParameters([], task, [parameters(1), parameters(2)]);
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
w = parameters(6);

%% Start loop over trials

p_initial = .5;
p_bias = .5; % Conservative bias towards a probability of .5 (i.e., neutral criterion)
z_model = zeros(1,NumTrials);
PChatA = z_model;
log_P = z_model;
beta = z_model;
p_estimate = zeros(NumTrials,1);
v_estimate = zeros(NumTrials,1);
p_conservative = p_estimate;
p_estimate(1,:) = p_initial;
p_conservative(1,:) = w*p_estimate(1,:) + (1-w)*p_bias;
beta(1) = p_conservative(1,:)/(1-p_conservative(1,:));
z_model(:,1) = sigma^2 * log(beta(1)) / diff(mu);
if task == 1
    log_P(:,1) = -0.5*log(2*pi*sigma_criterion) - 0.5*((resp_obs(1)-z_model(:,1))./sigma_criterion).^2;
else
    PChatA(:,1) = 1 - 0.5*erfc( -(S(1) - z_model(1)) / (sqrt(2)*sigma_ellipse) );   % Faster implementation
    log_P(:,1) = log(PChatA(:,1)).*(resp_obs(1)==1) + log(1-PChatA(:,1)).*(resp_obs(1)~=1);
end
Cprev = C(1);
    
% Step 1: Setup prior grid over r (reward probability - category probability in our task), 
% v (volatility), k (volatility's rate of change)

r = linspace(0.01, 0.99, 30);
v = linspace(log(0.01), log(0.99), 30);
k = linspace(log(5e-4), log(20), 30);

prior_r = 1/(r(end)-r(1));
prior_v = 1/(v(end)-v(1));
prior_k = 1/(k(end)-k(1));

prior = ones(numel(r), numel(v), numel(k)).*(prior_r*prior_v*prior_k);

delta_r = r(2)-r(1);
delta_v = v(2)-v(1);
delta_k = k(2)-k(1);

delta_grid = delta_r*delta_v*delta_k;

% Step 2: Define r and v transition probability matrices

% r-transition matrix
[r_1, r_i, v_1] = meshgrid(r, r, v);
r_trans = r_trans_func(r_1, r_i, v_1);
% normalize
temp = r_trans; for i_dim=1:3; temp = qtrapz(temp); end
z_norm = temp*delta_r*delta_v*delta_r;
r_trans = r_trans ./ z_norm;

% v-transition matrix
[k_i, v_i, v_1] = meshgrid(k, v, v);
v_trans = v_trans_func(k_i, v_i, v_1);
% normalize
temp = v_trans; for i_dim=1:3; temp = qtrapz(temp); end
z_norm = temp*delta_v*delta_v*delta_k;
v_trans = v_trans ./ z_norm;

for t = 2:NumTrials
    
    % Step 3: Multiply posterior probability from previous trial by p(v_current
    % | v_previous, k) and marginalize over v
    v_leaked = einsum('ijk,jkm->ikm', prior, v_trans);

    % Step 4: Multiple result by p(r_current | r_previous, v) and
    % marginalize over r
    r_leaked = einsum('iml,ikm->lmk', r_trans, v_leaked);
    
    % Step 5: Multiply by p(y_current | r_current)
    post_prob = zeros(numel(r), numel(v), numel(k));
    
    for idx_r1 = 1:numel(r)
        r_current = r(idx_r1);
        p_outcome = r_current^(Cprev)*(1-r_current)^(1-Cprev);
        for idx_v1 = 1:numel(v)
            for idx_k = 1:numel(k)
                post_prob(idx_r1, idx_v1, idx_k) = p_outcome*r_leaked(idx_r1, idx_v1, idx_k);
            end
        end
    end
    
    % normalize
    temp = post_prob; for i_dim=1:3; temp = qtrapz(temp); end
    z_norm = temp*delta_grid;
    post_prob = post_prob ./ z_norm;

    % Step 6: Marginalize v_current and k to get marginal distribution for
    % r_current
    prob_r1 = delta_v*delta_k*squeeze(qtrapz(qtrapz(post_prob, 2), 3));

    % Step 7a: Compute estimate of r_current
    p_estimate(t,:) = qtrapz(r.*prob_r1')*delta_r;

    % Step 7b: Compute estimate of v_current
    marg_v1 = delta_r*delta_k*squeeze(qtrapz(qtrapz(post_prob, 1), 3));
    v_estimate(t,:) = qtrapz(v.*marg_v1)*delta_v;

    % Step 8 (optional): Add conservatism bias
    p_conservative(t,:) = w*p_estimate(t,:) + (1-w)*p_bias;

    % Step 9: Compute the decision criterion
    beta(t) = p_conservative(t,:)/(1-p_conservative(t,:));
    z_model(t) = sigma^2 * log(beta(t)) / diff(mu);

    % Step 10: Compute log likelihood
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
    % update most recently observed data point and set prior to post
    Cprev = C(t);
    prior = post_prob;
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

%--------------------------------------------------------------------------
function [v_trans] = v_trans_func(k_i, v_i, v_1)
% V_TRANS_FUNC Volatility transition probability
    norm_constant = 1/sqrt(2*pi*exp(2*k_i));
    v_trans = norm_constant .* exp(-.5 * (v_1 - v_i).^2 ./ exp(2*k_i));
end

function [r_trans] = r_trans_func(r_1, r_i, v_1)
% R_TRANS_FUNC Category probability transition probability
    r_unif = 1/size(r_1,1);
    r_trans = (1 - exp(v_1)) .* (1 - (r_1 == r_i)) + exp(v_1).*r_unif;
end

end