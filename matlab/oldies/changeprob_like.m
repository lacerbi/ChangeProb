function [post,P,p_true,C,last,rmse] = changeprob_like(theta,data,method)
%CHANGEPROB_LIKE

if nargin < 1; theta = []; end
if nargin < 2 || isempty(data); data = 0; end
if nargin < 3 || isempty(method); method = 0; end

if isnumeric(data); rng(data); data = []; end
do_plots = 0;

%% Experiment constants
Ntrials = 800;
prior_rl = [80,120];    % Run length ~ Uniform[prior_rl(1),prior_rl(2)]
T = diff(prior_rl)+1;   % Number of elements
Nprobs = 5;
p_vec = linspace(0.2,0.8,Nprobs);

%% Simulated parameters
if isempty(theta)
    sigmaell = 10;
end

if ~isempty(data)
    %% Read session parameters
    col = data.SessionOrder(1);     % Column of overt criterion task
    
    sigma_gauss = data.StdDev;
    sigmatilde = sqrt(sigma_gauss^2 + sigmaell^2);
    mu = [data.GreenMean(col),data.RedMean(col)];
    mu = mu - mean(mu);
    C = (data.Category(:,col) == 2);    % Category A responses
    % C(C==0) = 2;                        % Category B
    p_true = data.pA(:,col)';
    X = sigmatilde*randn(1,Ntrials);    % Vector of noisy measurements

else
    %% Create fake dataset
    sigma_gauss = 10;
    sigmatilde = sqrt(sigma_gauss^2 + sigmaell^2);
    mu = 10;

    C = []; p_true = []; p = 0;
    while numel(C) < Ntrials
        runlength = randi(prior_rl);
        pold = p;
        while p == pold; p = p_vec(randi(Nprobs)); end
        p_true = [p_true, p*ones(1,runlength)];
        C = [C, 2 - (rand(1,runlength) < p)];          % Vector of true categories
    end
    C = C(1:Ntrials);
    X = sigmatilde*randn(1,Ntrials);    % Vector of noisy measurements
    p_true = p_true(1:Ntrials);
end

%% Initialize inference

tic

% priorp = bsxfun(@rdivide, priorp, sum(priorp,3));
%priorp = bsxfun(@rdivide, priorp, sum(priorp(:)));

% Hazard function (defined only where nonzero)
H = ones(T,1)/T;
H = H./flipud(cumsum(flipud(H)));

switch method
    case 0  % Old algorithm (something slightly off)
        % Posterior over run lengths (from 0 to PRIOR_RL(2)-1)
        post = zeros(prior_rl(2),Nprobs,Nprobs);
        post(1,:,1) = 1;  % Change in the first trial
        post = post ./ sum(post(:));
        priorp(1,:,:) = (ones(Nprobs,Nprobs) - eye(Nprobs)) / (Nprobs-1);
    
    case 1  % Bayesian OCPD extended
        % Posterior over run lengths (from 0 to PRIOR_RL(2)-1)
        post = zeros(prior_rl(2),Nprobs);
        post(1,:) = 1;  % Change in the first trial
        post = post ./ sum(post(:));
        tab = NaN(prior_rl(2),2); tab(1,:) = 0;
        priorp = (ones(Nprobs,Nprobs) - eye(Nprobs)); % /(Nprobs-1);
        
    case 2  % Bayesian OCPD fast
        % Posterior over run lengths (from 0 to PRIOR_RL(2)-1)
        post = zeros(prior_rl(2),Nprobs);
        post(1,:) = 1;  % Change in the first trial
        post = post ./ sum(post(:));
        tab = zeros(prior_rl(2),1,Nprobs); tab(1,1,:) = 1;
        priorp = (ones(Nprobs,Nprobs) - eye(Nprobs)); % /(Nprobs-1);
        p_vec3(1,1,:) = p_vec;
        priorp3(1,:,:) = priorp;
end

P = zeros(Ntrials,Nprobs);
last = zeros(Ntrials,prior_rl(2));

for t = 1:Ntrials
    %t
    switch method
        case 0
            post = bayesianOCPDupdate(X(t),C(t),post,[],priorp,H,p_vec,prior_rl,mu,sigmatilde,t);
            pp = sum(sum(post,3),1);
            tt = sum(sum(post,2),3);
        case 1
            [post,tab,temp] = bayesianOCPDupdate_explicit(X(t),C(t),post,tab,priorp,H,p_vec,prior_rl,mu,sigmatilde,t);
            pp = nansum(temp,1);
            tt = nansum(post,2);        
        case 2
            [post,tab,temp] = bayesianOCPDupdate_fast(X(t),C(t),post,tab,priorp3,H,p_vec3,prior_rl,mu,sigmatilde,t);
            pp = nansum(temp,1);
            tt = nansum(post,2);        
    end
    P(t,:) = pp/sum(pp);
    last(t,:) = tt/sum(tt);
end

toc

rmse = sqrt(mean((sum(bsxfun(@times,p_vec,P),2) - p_true').^2));

%% Plotting

if do_plots
    h = plotify([1 1, 2:Nprobs+1]','Margins',[0.05 0.05, 0.15 0.05],'Title',['Posterior over category probability for the BOCPD ideal observer']);
    axes(h(1));
    plot(1:Ntrials,p_true,'k','LineWidth',2);
    ylim([0,1]);
    box off;
    set(gca,'XTickLabel',[]);
    ylabel('Pr(A)');    
    text(Ntrials*0.01,0.9,['RMSE = ' num2str(rmse,'%.3g')]);
    
    change_vec = [0,find(diff(p_true) ~= 0),Ntrials];
    
    for i = 1:Nprobs
        axes(h(i+1));
        plot(1:Ntrials,P(:,i),'k'); hold on;
        box off;
        if i < Nprobs; set(gca,'XTickLabel',[]); end
        text(Ntrials*0.01,1.1,['Pr(A) = ' num2str(p_vec(i))])
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
    stdfig();
    xlabel('Trial number');
end

end

%--------------------------------------------------------------------------
function [post,tab] = bayesianOCPDupdate(x,C,post,tab,priorp,H,p_vec,prior_rl,mu,sigmatilde,t)
%BAYESIANCPDUPDATE Bayesian online changepoint detection update

    pxc1 = normpdf(x, mu, sigmatilde);
    pxc2 = normpdf(x, -mu, sigmatilde);

    if C == 1
        pc(1,:,1) = p_vec;
    else
        pc(1,:,1) = 1 - p_vec;
    end
    post = bsxfun(@times, post, pc);
        
    % Slice out trials where change is possible
    slice = post(prior_rl(1):prior_rl(2),:,:);
    
    % Multiply by Hazard function
    % slice = bsxfun(@times, slice, H);

    temp0 = sum(slice,3);
    temp2(:,1,:) = temp0;
    % temp = bsxfun(@times, temp2, priorp);
    % temp = temp2;
    temp = bsxfun(@times, temp2, H);

    currenttrial = sum(sum(temp,1),2);
    
    priorp3(1,:,:) = priorp;

    %if mod(t,100) == 0
    %    t
    %end
    
    slice = bsxfun(@times, slice, 1-H);
    post(prior_rl(1):prior_rl(2),:,:) = slice; 
    
    post = circshift(post,1,1);
    post(1,:,:) = repmat(currenttrial,[1,5,1]) .*priorp3;

    % Calculate evidence
    Px = sum(post(:));
    
    % Normalize posterior
    post = post ./ Px;

    % joint1 = joint1 - logsumexp(joint1(:));

    % joint1

end

%--------------------------------------------------------------------------
function [post,tab,postalpha] = bayesianOCPDupdate_explicit(x,C,post,tab,priorp,H,p_vec,prior_rl,mu,sigmatilde,t)
%BAYESIANCPDUPDATE Bayesian online changepoint detection update

    pxc1 = normpdf(x, mu, sigmatilde);
    pxc2 = normpdf(x, -mu, sigmatilde);

    idxrange = prior_rl(1):prior_rl(2);
    
    % 3. Evaluate predictive probability for C
    [pc,postpi] = predC(C,p_vec,tab,priorp);
        
    % (4-5) Multiply growth and changepoint probabilities by predictive prob.
    post = bsxfun(@times, post, pc);
    
    % 5. Calculate Changepoint Probabilities
    
    % Slice out trials where change is possible
    slice = post(idxrange,:);
    
    % Multiply by Hazard function
    % slice = bsxfun(@times, slice, H);
    
    % Compute posterior over beta_t
    p_vec3(1,:) = p_vec;
    pb = bsxfun(@power, p_vec3, tab(idxrange,1)) .*  bsxfun(@power, 1 - p_vec3, tab(idxrange,2));
    pb = bsxfun(@rdivide, pb, sum(pb,2));
        
    if mod(t,100) == 0
        t
    end
    
    prior3(1,:,:) = priorp;
    
    pt(:,1,:) = pb;
    temp = bsxfun(@times, pt, prior3);
    temp = bsxfun(@rdivide, pt, sum(pt,3));
    % temp = bsxfun(@times, pb, 1-bsxfun(@rdivide,slice_pre,nansum(slice_pre,2)));
    % temp = bsxfun(@rdivide, temp, sum(sum(temp,2),3));
    
    % currenttrial = nansum(bsxfun(@times, nansum(bsxfun(@times, slice, H),2), temp),1);
    currenttrial = nansum(sum(bsxfun(@times, bsxfun(@times, slice, H), temp),2),1);

    % currenttrial = nansum(bsxfun(@times,bsxfun(@times,sum(slice,2),H),temp),1);
    % currenttrial = nansum(bsxfun(@times,slice.*temp,H),1);
    % currenttrial = sum(sum(bsxfun(@times,slice,H),1),2);
    % currenttrial(isnan(currenttrial)) = -Inf;
    
    % Calculate growth probabilities    
    post(idxrange,:) = bsxfun(@times, slice, 1-H); 
    
    post = circshift(post,1,1);
    post(1,:) = currenttrial;

    post(isnan(post)) = 0;
    
    % Calculate evidence
    Px = sum(post(:));
    
    % Normalize posterior
    post = post ./ Px;

    % Update sufficient statistics
    tab = circshift(tab,1,1);
    tab(1,:) = 0;
    if C == 1
        tab(:,1) = tab(:,1) + 1;
    else
        tab(:,2) = tab(:,2) + 1;
    end
    
    % joint1 = joint1 - logsumexp(joint1(:));

    % joint1
    
    % Compute posterior over beta_t
    p_vec3(1,:) = p_vec;
    pb = bsxfun(@power, p_vec3, tab(:,1)) .*  bsxfun(@power, 1 - p_vec3, tab(:,2));
    pb = bsxfun(@rdivide, pb, nansum(pb,2));
    pt2(:,1,:) = pb;
    prior3(1,:,:) = priorp;    
    postalpha(:,:) = nansum(bsxfun(@times,bsxfun(@times, pt2, prior3), post),2);
    % temp = bsxfun(@times, pb, 1-bsxfun(@rdivide,slice_pre,nansum(slice_pre,2)));
    % temp = bsxfun(@rdivide, temp, sum(sum(temp,2),3));
    postalpha(isnan(postalpha)) = 0;


end

%--------------------------------------------------------------------------
function [post,tab,postalpha] = bayesianOCPDupdate_fast(x,C,post,tab,priorp,H,p_vec,prior_rl,mu,sigmatilde,t)
%BAYESIANCPDUPDATE Bayesian online changepoint detection update

    pxc1 = normpdf(x, mu, sigmatilde);
    pxc2 = normpdf(x, -mu, sigmatilde);

    idxrange = prior_rl(1):prior_rl(2);
    
    % 3. Evaluate predictive probability for C
    pi_post = bsxfun(@times, tab, priorp);               % Posterior over pi_i
    pi_post = bsxfun(@rdivide, pi_post, sum(pi_post,3));
    pc(:,:) = sum(bsxfun(@times, pi_post, p_vec),3);     % Predictive probability
    if C ~= 1; pc = 1-pc; end
        
    % (4-5) Multiply growth and changepoint probabilities by predictive prob.
    post = bsxfun(@times, post, pc);
    
    % 5. Calculate Changepoint Probabilities
    
    % Slice out trials where change is possible
    slice = post(idxrange,:);

    if mod(t,100) == 0
        t
    end
    
    % Posterior over pi_{t-1} that will become posterior over beta_{t}
    currenttrial = nansum(sum(bsxfun(@times, ...
        bsxfun(@times, slice, H), ...
        tab(idxrange,1,:)),2),1);
    
    % Calculate growth probabilities    
    post(idxrange,:) = bsxfun(@times, slice, 1-H); 
    
    post = circshift(post,1,1);
    post(1,:) = currenttrial;
    post(isnan(post)) = 0;
    
    % Calculate evidence
    Px = sum(post(:));
    
    % Normalize posterior
    post = post ./ Px;

    % Update sufficient statistics
    tab = circshift(tab,1,1);
    tab(1,1,:) = 1;
    if C == 1
        tab = bsxfun(@times, tab, p_vec);
    else
        tab = bsxfun(@times, tab, 1 - p_vec);
    end
    tab = bsxfun(@rdivide, tab, sum(tab,3));
    
    % Compute posterior over pi_t
    postalpha(:,:) = sum(bsxfun(@times,bsxfun(@times, tab, priorp), post),2);
    % postalpha(isnan(postalpha)) = 0;


end


function [p,postpi] = predC(C,p_vec,tab,priorp)
%PREDC Predictive probability for C

% Compute posterior over pi_i
p_vec3(1,:) = p_vec;
postpi = bsxfun(@power, p_vec3, tab(:,1)) .* bsxfun(@power, 1 - p_vec3, tab(:,2));
postpi = bsxfun(@rdivide, postpi, nansum(postpi,2));

pt(:,1,:) = postpi;
priorp3(1,:,:) = priorp;
p_vec4(1,1,:) = p_vec;

pt = bsxfun(@times, pt, priorp3);
pt = bsxfun(@rdivide, pt, sum(pt,3));

% Compute predictive probability
% p = bsxfun(@times, postpi, p_vec3) * priorp
p(:,:) = sum(bsxfun(@times, pt, p_vec4),3);

if C ~= 1; p = 1-p; end

end
