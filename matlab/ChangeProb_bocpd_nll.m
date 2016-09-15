function [nLL,post,P,p_true,C,last,rmse] = ChangeProb_bocpd_nll(parameters,data,task,options)
%CHANGEPROB_BOCPD_NLL Bayesian online changepoint detection observer.
% (Documentation to be written.)
%
% Author:   Luigi Acerbi
% Email:    luigi.acerbi@gmail.com
% Date:     09/15/2016

% Parameter vector:
% #1 is SIGMA_ELLIPSE, #2 is SIGMA_CRITERION, #3 is LAPSE rate
if nargin < 1; parameters = []; end

% Data struct or random seed for fake data generation
if nargin < 2 || isempty(data); data = 0; end
if isnumeric(data); rng(data); data = []; end

% Task (1 overt, 2 covert)
if nargin < 3 || isempty(task); task = 1; end
if task ~= 1 && task ~= 2; error('TASK can only be 1 (overt) or 2 (covert).'); end

% Additional options to set experimental constants
if nargin < 4; options = []; end

do_plot = nargout == 0; % If no outputs, make a plot

%% Experiment constants

% Run length ~ Uniform[prior_rl(1),prior_rl(2)]
if isfield(options,'prior_rl') && ~isempty(options.prior_rl)
    prior_rl = options.prior_rl;
else
    prior_rl = [80,120];        % Default runlengths
end
L = diff(prior_rl)+1;           % Number of elements

% Probability states
if isfield(options,'p_vec') && ~isempty(options.p_vec)
    p_vec = options.p_vec;
else
    p_vec = linspace(0.2,0.8,5);    % Default states
end
Nprobs = numel(p_vec);              % # states

%% Observer model parameters
switch numel(parameters)
    case 0  % Empty parameter vector
        sigma_ellipse = []; % Use SIGMA_ELLIPSE from calibration data
        sigma_criterion = [];
        lambda = 0;
    case 1  % No lapse rate
        sigma_ellipse = parameters(1);
        sigma_criterion = [];
        lambda = 0;
    case 2
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = 0;
    case 3
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = parameters(3);
end

if ~isempty(data)
    %% Read session parameters
    NumTrials = data.NumTrials;
    
    col = data.SessionOrder(1);     % Column of overt criterion task
    
    sigma_s = data.StdDev;
    if isempty(sigma_ellipse) || ~isfinite(sigma_ellipse); sigma_ellipse = data.EllipseNoise; end
    if isempty(sigma_criterion) || ~isfinite(sigma_criterion); sigma_criterion = sigma_ellipse; end
    sigma = sqrt(sigma_s^2 + sigma_ellipse^2);  % Category + sensory noise
    mu = [data.GreenMean(col),data.RedMean(col)];
    mu_bar = mean(mu);
    mu = mu - mu_bar;
    C = (data.Category(:,col) == 2);        % Category A
    p_true = data.pA(:,col);
    X = sigma*randn(1,NumTrials);    % Vector of noisy measurements
    
    beta_resp = data.Criterion(:,col) - mu_bar;   % Reported criterion

else
    %% Create fake dataset
    NumTrials = 80;
    sigma_s = 10;
    if isempty(sigma_ellipse) || ~isfinite(sigma_ellipse); sigma_ellipse = 5.68; end
    if isempty(sigma_criterion) || ~isfinite(sigma_criterion); sigma_criterion = sigma_ellipse; end
    sigma = sqrt(sigma_s^2 + sigma_ellipse^2);
    mu = [-10,10];

    C = []; p_true = []; p = 0;
    while numel(C) < NumTrials
        runlength = randi(prior_rl);
        pold = p;
        while p == pold; p = p_vec(randi(Nprobs)); end
        p_true = [p_true; p*ones(runlength,1)];
        C = [C; 2 - (rand(runlength,1) < p)];          % Vector of true categories
    end
    C = C(1:NumTrials);
    p_true = p_true(1:NumTrials);
    
    % Generate stimuli based on category labels
    S = mu(C)' + sigma_s*randn(NumTrials,1);
    
    % Generate noisy measurements
    X = S + sigma_ellipse*randn(NumTrials,1);
    
    beta_resp = zeros(NumTrials,1);
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

% Current posterior probability over state
pi_post = ones(1,Nprobs)/Nprobs;

%% Begin loop over trials
P = zeros(NumTrials+1,Nprobs);
P(1,:) = ones(1,Nprobs)/Nprobs;
last = zeros(NumTrials,size(post,1));

for t = 1:NumTrials
    %t
    [post,Psi,pi_post] = bayesianOCPDupdate(X(t),C(t),post,Psi,Tmat,H,p_vec3,mu,sigma,t);
    tt = nansum(post,2);        

    % If the posterior is about the next trial, this should be t+1?
    P(t,:) = pi_post;
    %P(t+1,:) = pi_post;
    last(t,:) = tt/sum(tt);
end

% Compute predicted criterion for the overt task
Gamma_t = sum(bsxfun(@times, P, p_vec(:)'),2) ./ sum(bsxfun(@times, P, 1-p_vec(:)'), 2);
beta_opt = sigma^2 * log(Gamma_t) / diff(mu);

% Log probability of overt task responses
log_Pbeta = -0.5*log(2*pi*sigma_criterion) - 0.5*((beta_resp-beta_opt(1:NumTrials))./sigma_criterion).^2;
if lambda > 0
    log_Pbeta = log(lambda/360 + (1-lambda)*exp(log_Pbeta));    
end

% Sum log likelihood
nLL = -nansum(log_Pbeta);

meanP = sum(bsxfun(@times,p_vec,P(1:NumTrials,:)),2);

rmse = sqrt(mean((meanP - p_true).^2));

%rmse = sqrt(mean((sum(bsxfun(@times,p_vec,P(1:end-1,:)),2) - p_true(2:end)').^2));
%[rmse,rmse2]

%% Plotting
if do_plot
    h = plotify([1 1, 2:Nprobs+1, (Nprobs+2)*[1 1 1]]','Margins',[0.05 0.05, 0.15 0.05],'Title',['Posterior over category probability for the BOCPD ideal observer']);
    axes(h(1));
    plot(1:NumTrials,p_true,'k','LineWidth',2); hold on;
    plot(1:NumTrials,meanP,'r','LineWidth',2);
    ylim([0,1]);
    box off;
    set(gca,'XTickLabel',[]);
    ylabel('Pr(A)');    
    text(NumTrials*0.01,0.9,['RMSE = ' num2str(rmse,'%.3g')]);
    
    change_vec = [0;find(diff(p_true) ~= 0);NumTrials];
    
    for i = 1:Nprobs
        axes(h(i+1));
        plot(1:NumTrials,P(1:NumTrials,i),'k'); hold on;
        box off;
        if i < Nprobs; set(gca,'XTickLabel',[]); end
        text(NumTrials*0.01,1.1,['Pr(A) = ' num2str(p_vec(i))])
        for j = 1:numel(change_vec)-1
            plot(change_vec(j)*[1 1],[0 1],':k','LineWidth',1);
            if abs(p_true(change_vec(j)+1) - p_vec(i)) < 1e-8
            % if p_true(change_vec(j)+1) == p_vec(i)
                xx = [change_vec(j) change_vec(j+1) change_vec(j+1) change_vec(j)];
                yy = [0 0, 1 1];
                fill(xx,yy,'c','EdgeColor','none','FaceAlpha',0.5);                
            end
        end
        
    end
    
    axes(h(Nprobs+2));
    
    % Smoothen responded criterion for visualization
    beta_smooth = smooth(beta_resp,11);
    
    scatter(1:NumTrials,beta_resp,'k.'); hold on;
    plot(1:NumTrials,beta_smooth,'k','LineWidth',2); hold on;
    plot(1:NumTrials,beta_opt(1:NumTrials),'r','LineWidth',2);
    ylabel('Criterion');
    ylim([-50 120]);
    set(gca,'YTick',[-50,0,50]);
    h_leg = legend('Reported criterion','Moving average','Model');
    set(h_leg,'Box','off','Location','NorthEast');
    % set(get(h_leg,'Children'),'Visible','on');
    
    stdfig();
    xlabel('Trial number');
end

end

%--------------------------------------------------------------------------
function [post,Psi,pi_post] = bayesianOCPDupdate(x,C,post,Psi,Tmat,H,p_vec,mu,sigma,t)
%BAYESIANCPDUPDATE Bayesian online changepoint detection update

    % This FLAG sets how the algorithm performs some computations.
    % The correct way should be FLAG=0 but the algorithm performs better 
    % with FLAG=1, which suggests something might be off.
    flag = 0;

    %pxc1 = normpdf(x, -mu, sigma);
    %pxc2 = normpdf(x, mu, sigma);

    %if mod(t,100) == 0
    %    t
    %end
    L = size(H,1);  % Width of prior over run lengths
    % Slice with nonzero hazard function (probability of change)
    idxrange = (size(post,1)-L+1:size(post,1));
    
    % 3. Evaluate predictive probability
    pi_post = bsxfun(@times, Psi, Tmat);                    % Posterior over pi_i
    pi_post = bsxfun(@rdivide, pi_post, sum(pi_post,3));
    predC(:,:) = sum(bsxfun(@times, pi_post, p_vec),3);     % Predictive probability
    if C ~= 1; predC = 1-predC; end
        
    % 3b. Multiply posterior by predictive probability
    post = bsxfun(@times, post, predC);
    
    % 4. Calculate Changepoint Probabilities
    
    % Slice out trials where change is possible
    slice = post(idxrange,:);
        
    % Posterior over pi_{t-1} that will become posterior over xi_t
    if flag
        currenttrial = nansum(sum(bsxfun(@times, ...
            bsxfun(@times, slice, H), ...
            Psi(idxrange,1,:)),2),1);
    else
        currenttrial = nansum(sum(bsxfun(@times,bsxfun(@times, ...
            bsxfun(@times, slice, H), ...
            Psi(idxrange,1,:)),Tmat),2),1);
    end
    
    % 5. Calculate growth probabilities    
    post(idxrange,:) = bsxfun(@times, slice, 1-H); 
    
    post = circshift(post,1,1);
    post(1,:) = currenttrial;
    post(isnan(post)) = 0;
    
    % 6. Calculate normalization
    Z = sum(post(:));
    
    % 7. Normalize posterior
    post = post ./ Z;

    % 8. Update sufficient statistics
    Psi = circshift(Psi,1,1);
    if flag; Psi(1,1,:) = 1; end
    if C == 1
        Psi = bsxfun(@times, Psi, p_vec);
    else
        Psi = bsxfun(@times, Psi, 1 - p_vec);
    end
    if ~flag; Psi(1,1,:) = 1; end
    Psi = bsxfun(@rdivide, Psi, sum(Psi,3));
    
    % Compute marginalized posterior over pi_t
    pi_post = nansum(sum(bsxfun(@times,bsxfun(@times, Psi, Tmat), post),2),1);
    pi_post = pi_post / sum(pi_post);

end

%--------------------------------------------------------------------------
function stdfig(fig)
%STDFIG Standardized figure.
if nargin < 1 || isempty(fig); fig = gcf; end
screensize = get(groot,'ScreenSize');
c = fig.Children;
for i = 1:numel(c)
   if ~isa(c(i),'matlab.graphics.axis.Axes'); continue; end 
   set(c(i),'Box','off','TickDir','out');    
end
set(fig,'Color','w');
set(fig,'Position',screensize);
end