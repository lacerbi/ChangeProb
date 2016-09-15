function [joint1,P,p_true,C,last,rmse] = changeprob_like(K)
%CHANGEPROB_LIKE 

if nargin < 1 || isempty(K); K = 1; end

rng(0);
Ntrials = 800;
do_plots = 1;

%% Model parameters and constants

sigma = 1;
sigmas = 1;
sigmatilde = sqrt(sigma^2 + sigmas^2);
mu = 1;

changerange = [80,120]; % Change happens in a range every # trials
T = diff(changerange)+1;
Nprobs = 5;
p_vec = linspace(0.2,0.8,Nprobs);

%% Create fake dataset
x = sigmatilde*randn(1,Ntrials);    % Vector of noisy measurements
C = []; p_true = []; p = 0;
while numel(C) < Ntrials
    runlength = randi(changerange);
    pold = p;
    while p == pold; p = p_vec(randi(Nprobs)); end
    p_true = [p_true, p*ones(1,runlength)];
    C = [C, 2 - (rand(1,runlength) < p)];          % Vector of true categories
end
C = C(1:Ntrials);
p_true = p_true(1:Ntrials);

%% Initialize inference

tic

switch K
    case 1
        priorp(1,:,:) = ones(Nprobs,Nprobs) - eye(Nprobs);
        priorp = bsxfun(@rdivide, priorp, sum(priorp,3));
        joint1 = zeros(changerange(2),Nprobs,Nprobs);
        joint1(1,:,1) = 1/Nprobs;  % Change in the first trial
    case 2
        priorp(1,1,:,:,:) = ones(Nprobs,Nprobs,Nprobs);
        for i = 1:Nprobs
            for j = 1:Nprobs
                for k = 1:Nprobs
                    if i == j || j == k; priorp(1,1,i,j,k) = 0; end
                end
            end
        end
        priorp = bsxfun(@rdivide, priorp, sum(priorp,5));
        joint1 = zeros(changerange(2),T,Nprobs,Nprobs,Nprobs);
        joint1(1,1,:,:,1) = ones(Nprobs,Nprobs) - eye(Nprobs);  % Change in the first trial
end

joint1 = joint1 ./ sum(joint1(:));
joint1 = log(joint1);

for t = 1:Ntrials
    %t
    probC_nexttrial(K);
    if K == 1
        pp = logsumexp(logsumexp(joint1,3),1);
        pp = exp(pp-max(pp));
        P(t,:) = pp/sum(pp);
        tt = logsumexp(logsumexp(joint1,2),3);
        last(t,:) = exp(tt)/exp(logsumexp(tt));
    else
        pp = logsumexp(logsumexp(logsumexp(logsumexp(joint1,5),4),2),1);
        pp = exp(pp-max(pp));
        P(t,:) = pp/sum(pp);
        tt = logsumexp(logsumexp(logsumexp(logsumexp(joint1,2),3),4),5);
        last(t,:) = exp(tt)/exp(logsumexp(tt));        
    end
end

toc

rmse = sqrt(mean((sum(bsxfun(@times,p_vec,P),2) - p_true').^2));

%% Plotting

if do_plots
    h = plotify([1 1, 2:Nprobs+1]','Margins',[0.05 0.05, 0.15 0.05],'Title',['Posterior over category probability for the truncated ideal observer (K = ' num2str(K) ')']);
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
            if p_true(change_vec(j)+1) == p_vec(i)
                xx = [change_vec(j) change_vec(j+1) change_vec(j+1) change_vec(j)];
                yy = [0 0, 1 1];
                fill(xx,yy,'c','EdgeColor','none','FaceAlpha',0.5);                
            end
        end
        
    end
    stdfig();
    xlabel('Trial number');
end



    function probC_nexttrial(K)

        pxc1 = normpdf(x(t), mu, sigmatilde);
        pxc2 = normpdf(x(t), -mu, sigmatilde);
        
        % Effective prior width
        % Teff = min(1+t-changerange(1),T);
        Teff = T;
        
        if K == 1
            lpc1(1,:,1) = log(p_vec);
            lpc2(1,:,1) = log(1-p_vec);
        
            % Joint slice from previous changepoints
            slice = joint1(changerange(1):changerange(2),:,:);

            temp0 = logsumexp(slice,3);

            temp2(:,1,:) = temp0;
            temp = bsxfun(@plus, temp2, log(priorp));
            temp(isnan(temp)) = -Inf;

            currenttrial = logsumexp(temp,1);
            if Teff > 0
                currenttrial = currenttrial - log(Teff);
            end
            currenttrial(isnan(currenttrial)) = -Inf;

            joint1 = circshift(joint1,1,1);
            joint1(1,:,:) = currenttrial;            
            
        else
            lpc1(1,1,:) = log(p_vec);
            lpc2(1,1,:) = log(1-p_vec);
            
            % Joint slice from previous changepoints
            slice = joint1(changerange(1):changerange(2),:,:,:,:);
            
            %if t == 400
            %    t
            %end
            
            temp0 = logsumexp(slice,5);
            temp2(:,:,1,:,:) = temp0;
            temp = bsxfun(@plus, temp2, log(priorp));
            temp(isnan(temp)) = -Inf;
            currenttrial = logsumexp(temp,2);
            if Teff > 0 
                currenttrial = currenttrial - log(Teff);
            end
            currenttrial(isnan(currenttrial)) = -Inf;

            joint1 = circshift(joint1,1,1);
            joint1(1,:,:,:,:) = currenttrial;            
            
        end
        
        if C(t) == 1
            joint1 = bsxfun(@plus, joint1, lpc1);
        else
            joint1 = bsxfun(@plus, joint1, lpc2);
        end

        % joint1 = joint1 - logsumexp(joint1(:));
        
        % joint1
        
    end

end