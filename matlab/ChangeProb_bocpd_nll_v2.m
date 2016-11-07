function [nLL, rmse, p_estimate, resp_model, post] = ChangeProb_bocpd_nll_v2(parameters, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task, prior_rl)
%CHANGEPROB_BOCPD_NLL Bayesian online changepoint detection observer.
% (Documentation to be written.)
%
% Author:   Luigi Acerbi
% Email:    luigi.acerbi@gmail.com
% Date:     Oct/1/2016

% Modified by Elyse Norton on Oct/27/2016

% Parameter vector:
% #1 is SIGMA_ELLIPSE, #2 is SIGMA_CRITERION, #3 is LAPSE, #4 is GAMMA,
% #5 is ALPHA, and #6 is W
if nargin < 1
    parameters = [];
    prior_rl = [];
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs] = changeprob_getSessionParameters();
    task = 1; 
end

% Data struct or random seed for fake data generation
if nargin < 2; error('You must specify the session parameters.'); end

% Task (1 overt, 2 covert)
if nargin < 9 || isempty(task); task = 1; end
if task ~= 1 && task ~= 2; error('TASK can only be 1 (overt-criterion) or 2 (covert-criterion).'); end

if nargin < 10 || isempty(prior_rl)
    prior_rl = [80,120];    % Default runlengths
end

%% Experiment constants

% Run length ~ Uniform[prior_rl(1),prior_rl(2)]
% if isfield(options,'prior_rl') && ~isempty(options.prior_rl)
%     prior_rl = options.prior_rl;
% else
%end
L = diff(prior_rl)+1;           % Number of elements

% Probability states
% if isfield(options,'p_vec') && ~isempty(options.p_vec)
%     p_vec = options.p_vec;
% else
    p_vec = linspace(0.2,0.8,5);    % Default states
%end
Nprobs = numel(p_vec);              % # states

%% Observer model parameters

switch numel(parameters)
    case 0  % Empty parameter vector
        sigma_criterion = sigma_ellipse;
        lambda = 0;
        gamma = Inf;
    case 1
        sigma_ellipse = parameters(1);
    case 2
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
    case 3
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = parameters(3);
    case 4
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = parameters(3);
        gamma = parameters(4);
    otherwise
        error('PARAMETERS should be a vector with up to four parameters.');
end

%% Initialize inference

% Hazard function (defined only where nonzero)
H = ones(L,1)/L;
H = H./flipud(cumsum(flipud(H)));

% Posterior over run lengths (from 0 to PRIOR_RL(2)-1)
post = zeros(prior_rl(2),Nprobs);
post(1,:) = 1;  % Change in the first trial
post = post ./ sum(post(:));

% Table of binomial count/probability
Psi = zeros(size(post,1),1,Nprobs); Psi(1,1,:) = 1/Nprobs;

% Transition matrix
Tmat(1,:,:) = (ones(Nprobs,Nprobs) - eye(Nprobs)) / (Nprobs-1);
% Tmat(1,:,:) = 1;
p_vec3(1,1,:) = p_vec;

% Auxiliary variables for covert-criterion task
if task == 2
    nx = 1001;  % Measurement grid size
    MAXSD = 5;  % When integrating a Gaussian go up to this distance
    X = bsxfun(@plus, S, sigma_ellipse*MAXSD*linspace(-1,1,nx));
    w = normpdf(linspace(-MAXSD,MAXSD,nx));
    w = w./qtrapz(w);
end

%% Begin loop over trials
P = zeros(NumTrials+1,Nprobs);
P(1,:) = ones(1,Nprobs)/Nprobs;
last = zeros(NumTrials,size(post,1));
PCx = [];

for t = 1:NumTrials
    %t
    if task == 2; Xt = X(t,:); else Xt = []; end    
    [post,Psi,pi_post,PCxA] = bayesianOCPDupdate(Xt,C(t),post,Psi,Tmat,H,p_vec3,mu,sigma,task);
    tt = nansum(post,2);

    % The predictive posterior is about the next trial
    P(t+1,:) = pi_post;
    last(t,:) = tt/sum(tt);
    
    % Record conditional posterior for covert-criterion task
    if task == 2
        if isempty(PCx); PCx = zeros(NumTrials,size(PCxA,2)); end
        PCx(t,:) = PCxA;
    end
end

%% Compute log likelihood
switch task
    case 1  % Overt-criterion task
        
        % Compute predicted criterion for the overt task
        Gamma_t = sum(bsxfun(@times, P, p_vec(:)'),2) ./ sum(bsxfun(@times, P, 1-p_vec(:)'), 2);
        z_opt = sigma^2 * log(Gamma_t) / diff(mu);

        % Log probability of overt task responses
        log_Pz = -0.5*log(2*pi*sigma_criterion) - 0.5*((resp_obs-z_opt(1:NumTrials))./sigma_criterion).^2;
        if lambda > 0
            log_Pz = log(lambda/360 + (1-lambda)*exp(log_Pz));    
        end

        % Sum negative log likelihood
        nLL = -nansum(log_Pz);
        resp_model = z_opt(1:NumTrials);
        
    case 2  % Covert-criterion task
        
        Pcx = 1./(1 + ((1-PCx)./PCx).^gamma);       % Softmax        
        PChatA = qtrapz(bsxfun(@times,Pcx,w),2);    % Marginalize over noise
        lambda = max(1e-4,lambda);                  % Minimum lapse to avoid numerical trouble
        PChatA = lambda/2 + (1-lambda)*PChatA;
        
        % Log probability of covert task responses
        log_PChat = log(PChatA).*(resp_obs == 1) + log(1-PChatA).*(resp_obs ~= 1);
        
        % Sum negative log likelihood
        nLL = -nansum(log_PChat);   
        resp_model = PChatA(1:NumTrials);
end

% RMSE between predictive posterior probability and true category probability
meanP = sum(bsxfun(@times,p_vec,P(1:NumTrials,:)),2);
p_estimate = meanP;
rmse = sqrt(mean((meanP - p_true).^2));

end

%--------------------------------------------------------------------------
function [post,Psi,pi_post,PCxA] = bayesianOCPDupdate(X,C,post,Psi,Tmat,H,p_vec,mu,sigma,task)
%BAYESIANCPDUPDATE Bayesian online changepoint detection update

    %if mod(t,100) == 0
    %    t
    %end
    L = size(H,1);  % Width of prior over run lengths
    % Slice with nonzero hazard function (probability of change)
    idxrange = (size(post,1)-L+1:size(post,1));

    % Posterior over pi_i before observing C_t
    pi_post = bsxfun(@times, Psi, Tmat);
    predCatA(:,:) = sum(bsxfun(@times, pi_post, p_vec),3);
    predCatB(:,:) = sum(bsxfun(@times, pi_post, 1 - p_vec),3);
    
    %----------------------------------------------------------------------
    % 2. Compute probability of response for covert-criterion task
    if task == 2
        pxCatA = exp(-0.5*((X-mu(1))./sigma).^2);
        pxCatB = exp(-0.5*((X-mu(2))./sigma).^2);

        PCxA = nansum(nansum(bsxfun(@times,predCatA,post),2),1).*pxCatA;
        PCxB = nansum(nansum(bsxfun(@times,predCatB,post),2),1).*pxCatB;    
        PCxA = PCxA./(PCxA + PCxB);
    else
        PCxA = [];
    end
    
    %----------------------------------------------------------------------
    % 3. Observe Ct (do nothing)
    
    %----------------------------------------------------------------------
    % 4a. Evaluate predictive probability
    if C == 1
        predCat = predCatA ./ (predCatA + predCatB);
    else
        predCat = predCatB ./ (predCatA + predCatB);         
    end
    
    % 4b. Multiply posterior by predictive probability
    post = bsxfun(@times, post, predCat);
        
    % 4c. Evaluate posterior probability over state (only for relevant range)
    if C == 1
        pi_postC = bsxfun(@times, pi_post(idxrange,:,:), p_vec);
    else
        pi_postC = bsxfun(@times, pi_post(idxrange,:,:), 1-p_vec);
    end
    pi_postC = bsxfun(@rdivide, pi_postC, sum(pi_postC,3));

    %----------------------------------------------------------------------    
    % 5a. Calculate unnormalized changepoint probabilities
    slice = post(idxrange,:);   % Slice out trials where change is possible
    currenttrial = nansum( ...
        bsxfun(@times, sum(bsxfun(@times, slice, pi_postC),2), H), ...
        1);
    
    % 5b. Calculate unnormalized growth probabilities    
    post(idxrange,:) = bsxfun(@times, slice, 1-H); 
    
    % Shift posterior by one step
    post = circshift(post,1,1);
    post(1,:) = currenttrial;
    post(isnan(post)) = 0;
    
    % 5c. Calculate normalization
    Z = sum(post(:));
    
    % 5d. Normalize posterior
    post = post ./ Z;
    
    %----------------------------------------------------------------------
    % 6a. Update sufficient statistics
    Psi = circshift(Psi,1,1);
    if C == 1
        Psi = bsxfun(@times, Psi, p_vec);
    else
        Psi = bsxfun(@times, Psi, 1 - p_vec);
    end
    Psi(1,1,:) = 1;
    Psi = bsxfun(@rdivide, Psi, sum(Psi,3));
    
    % 6b. Store predictive posterior over pi_t
    pi_post = nansum(sum(bsxfun(@times,bsxfun(@times, Psi, Tmat), post),2),1);
    pi_post = pi_post / sum(pi_post);

end
