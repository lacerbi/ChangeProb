function [vp,elbo,elbo_sd,exitflag] = vbmc_diagnostics(vps,outputs,beta_lcb)
%VBMC_DIAGNOSTICS Perform convergence diagnostics using outputs from VBMC.
%

if nargin < 3 || isempty(beta_lcb); beta_lcb = 3; end

if isstruct(outputs) && numel(outputs) == 1
    error('Wrong inputs.');
end

if isstruct(vps)
    for i = 1:numel(vps); temp{i} = vps(i); end
    vps = temp;
    clear temp;
end

if isstruct(outputs)
    for i = 1:numel(outputs); temp{i} = outputs(i); end
    outputs = temp;
    clear temp;    
end

Nfits = numel(outputs);
        
% Compute KL-divergence across all pairs of solutions
kl_mat = zeros(Nfits,Nfits);
for iRun = 1:Nfits
    for jRun = iRun+1:Nfits        
        kl = vbmc_kldiv(vps{iRun},vps{jRun});
        kl_mat(iRun,jRun) = kl(1);
        kl_mat(jRun,iRun) = kl(2);        
    end
end

kl_mat

for iFit = 1:Nfits
    elbo(iFit) = outputs{iFit}.elbo;
    elbo_sd(iFit) = outputs{iFit}.elbosd;
    exitflag(iFit) = strcmpi(outputs{iFit}.convergencestatus,'probable');
end

idx_ok = exitflag == 1;

if sum(idx_ok) == 0
    idx_ok = true(size(idx_ok)); 
    warning('No solution has converged, using potentially unstable solution.');
end

elcbo = elbo - beta_lcb*elbo_sd;
elcbo
elcbo(~idx_ok) = -Inf;
[~,best_idx] = max(elcbo);

vp = vps{best_idx};    
elbo = elbo(best_idx);
elbo_sd = elbo_sd(best_idx);
exitflag = exitflag(best_idx);

end