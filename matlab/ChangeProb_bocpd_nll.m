function [nLL,post,P,p_true,C,last,rmse] = ChangeProb_bocpd_nll(parameters,data,task,options)
%CHANGEPROB_BOCPD_NLL Bayesian online changepoint detection observer.
% (Documentation to be written.)
%
% Author:   Luigi Acerbi
% Email:    luigi.acerbi@gmail.com
% Date:     Oct/1/2016

% Parameter vector:
% #1 is SIGMA_ELLIPSE, #2 is SIGMA_CRITERION, #3 is LAPSE rate, #4 is GAMMA
if nargin < 1; parameters = []; end

% Data struct or random seed for fake data generation
if nargin < 2 || isempty(data); data = 0; end
if isnumeric(data); rng(data); data = []; end

% Task (1 overt, 2 covert)
if nargin < 3 || isempty(task); task = 1; end
if task ~= 1 && task ~= 2; error('TASK can only be 1 (overt-criterion) or 2 (covert-criterion).'); end

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
sigma_ellipse = [];     % Default: use SIGMA_ELLIPSE from calibration data
sigma_criterion = [];   % Default: use SIGMA_ELLIPSE from calibration data
lambda = 0;             % Default no lapse
gamma = Inf;            % Default deterministic covert-criterion choice

switch numel(parameters)
    case 0  % Empty parameter vector
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

if ~isempty(data)
    %% Read session parameters
    col = data.SessionOrder(task);     % Column of chosen task
    NumTrials = data.NumTrials;    
    
    % Noise parameters
    sigma_s = data.StdDev;
    if isempty(sigma_ellipse) || ~isfinite(sigma_ellipse); sigma_ellipse = data.EllipseNoise; end
    if isempty(sigma_criterion) || ~isfinite(sigma_criterion); sigma_criterion = sigma_ellipse; end
    
    % Category information 
    % (in the data Category B/Red is coded as 1, Category A/Green is coded as 2)
    C = (data.Category(:,col) == 2);        % Category A/Green
    p_true = data.pA(:,col);
    mu = [data.GreenMean(col),data.RedMean(col)];
    
    % Shift coordinate system to zero
    mu_bar = mean(mu);
    mu = mu - mu_bar;
    S = data.StimulusAngle(:,col) - mu_bar;
    
    % Get task-relevant responses
    switch task
        case 1  % Overt-criterion task    
            beta_resp = data.Criterion(:,col) - mu_bar;   % Reported criterion
        case 2  % Covert-criterion task
            Chat = data.Response(:,col) == 2;
    end

else
    %% Create fake dataset
    NumTrials = 800;
    sigma_s = 10;
    if isempty(sigma_ellipse) || ~isfinite(sigma_ellipse); sigma_ellipse = 5.68; end
    if isempty(sigma_criterion) || ~isfinite(sigma_criterion); sigma_criterion = sigma_ellipse; end
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
    
    beta_resp = zeros(NumTrials,1);
    Chat = NaN(NumTrials,1);
end

% Category + sensory noise
sigma = sqrt(sigma_s^2 + sigma_ellipse^2);

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
        beta_opt = sigma^2 * log(Gamma_t) / diff(mu);

        % Log probability of overt task responses
        log_Pbeta = -0.5*log(2*pi*sigma_criterion) - 0.5*((beta_resp-beta_opt(1:NumTrials))./sigma_criterion).^2;
        if lambda > 0
            log_Pbeta = log(lambda/360 + (1-lambda)*exp(log_Pbeta));    
        end

        % Sum negative log likelihood
        nLL = -nansum(log_Pbeta);
        
    case 2  % Covert-criterion task
        
        Pcx = 1./(1 + ((1-PCx)./PCx).^gamma);       % Softmax        
        PChatA = qtrapz(bsxfun(@times,Pcx,w),2);    % Marginalize over noise
        lambda = max(1e-4,lambda);                  % Minimum lapse to avoid numerical trouble
        PChatA = lambda/2 + (1-lambda)*PChatA;
        
        % Log probability of covert task responses
        log_PChat = log(PChatA).*(Chat == 1) + log(1-PChatA).*(Chat ~= 1);
        
        % Sum negative log likelihood
        nLL = -nansum(log_PChat);        
end

% RMSE between predictive posterior probability and true category probability
if do_plot || nargout > 6
    meanP = sum(bsxfun(@times,p_vec,P(1:NumTrials,:)),2);
    rmse = sqrt(mean((meanP - p_true).^2));
end

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
        
        for j = 1:numel(change_vec)-1
            plot(change_vec(j)*[1 1],[0 1],':k','LineWidth',1);
            if abs(p_true(change_vec(j)+1) - p_vec(i)) < 1e-8
            % if p_true(change_vec(j)+1) == p_vec(i)
                xx = [change_vec(j) change_vec(j+1) change_vec(j+1) change_vec(j)];
                yy = [0 0, 1 1];
                % fill(xx,yy,'c','EdgeColor','none','FaceAlpha',0.5);                
                fill(xx,yy,'c','EdgeColor','none'); hold on;
            end
        end
        
        plot(1:NumTrials,P(1:NumTrials,i),'k'); hold on;
        box off;
        if i < Nprobs; set(gca,'XTickLabel',[]); end
        text(NumTrials*0.01,1.1,['Pr(A) = ' num2str(p_vec(i))])
        
    end
    
    axes(h(Nprobs+2));
    
    switch task
        case 1
            % Smoothen responded criterion for visualization
            beta_smooth = smooth(beta_resp,11);
            
            if ~isempty(data)
                scatter(1:NumTrials,beta_resp,'k.'); hold on;
                plot(1:NumTrials,beta_smooth,'k','LineWidth',2); hold on;
            end
            plot(1:NumTrials,beta_opt(1:NumTrials),'r','LineWidth',2);
            ylabel('Criterion');
            ylim([-50 120]);
            set(gca,'YTick',[-50,0,50]);
            if ~isempty(data)
                h_leg = legend('Reported criterion','Moving average','Model');
            else
                h_leg = legend('Model');                
            end
            set(h_leg,'Box','off','Location','NorthEast');
            
        case 2
            % Smoothen responded categories for visualization
            Chat_smooth = smooth(double(Chat),11);
            if ~isempty(data)
                scatter(1:NumTrials,double(Chat),'k.'); hold on;
                plot(1:NumTrials,Chat_smooth,'k','LineWidth',2); hold on;
            end
            plot(1:NumTrials,smooth(PChatA(1:NumTrials),11),'r','LineWidth',2);
            ylabel('Probability response (C=A)');
            ylim([0 1.75]);
            set(gca,'YTick',[0 1]);
            if ~isempty(data)
                h_leg = legend('Reported category','Moving average (data)','Moving average (model)');
            else
                h_leg = legend('Moving average (model)');                
            end
            set(h_leg,'Box','off','Location','NorthEast');
    end
    
    stdfig();
    xlabel('Trial number');
end

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