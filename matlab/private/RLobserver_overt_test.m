%RLOBSERVER_OVERT_TEST
%   Test script for MEX-file RLobserver_overt.
%
%   Template MATLAB code generated on 05-Oct-2016 with MEXXER v0.2 
%   (https://github.com/lacerbi/mexxer).

TolErr = sqrt(eps);	% Maximum error tolerance per array element

% Define array sizes (choose reasonable values)
Ns = 1000;
Nt = 800;

% Randomly initialize input variables
% (or write here alternative initializations)
parameters = rand([1,2]);	%PARAMETERS: learning rate and adjustment noise.
sigma = 10*rand([1,1]);	%SIGMA: combined category and sensory noise.
dmu = 10*rand([1,1]);	%DMU: distance between category means.
X = 10*rand([Nt,Ns]);	%X: matrix of noisy measurements.
p_initial = rand([1,1]);	%P_INITIAL: initial probability.
z_resp = 10*rand([Nt,1]);	%Z_RESP: subject's criterion responses.
score = randi(2,[Nt,1])-1;	%SCORE: Trial feedback.

fprintf('=========================\n');
fprintf('Testing RLobserver_overt:\n');
fprintf('=========================\n');

% Call MATLAB and MEX functions
tic; [model_resp,log_Pz] = RLobserver_overt(parameters,sigma,dmu,X,p_initial,z_resp,score); t = toc;
tic; [model_resp_mex,log_Pz_mex] = RLobserver_overt_mex(parameters,sigma,dmu,X,p_initial,z_resp,score); t_mex = toc;

% Correctness check
model_resp_err = sum(abs(model_resp(:)-model_resp_mex(:)));
fprintf('Total error (model_resp): %g\n', model_resp_err);
if model_resp_err > TolErr*numel(model_resp);
	error('mexxer:tooLargeError','Correctness check failed. Error too large in model_resp.');
end
log_Pz_err = sum(abs(log_Pz(:)-log_Pz_mex(:)));
fprintf('Total error (log_Pz): %g\n', log_Pz_err);
if log_Pz_err > TolErr*numel(log_Pz);
	error('mexxer:tooLargeError','Correctness check failed. Error too large in log_Pz.');
end

% Runtime analysis
fprintf('Time for MATLAB code: %.3f s\n', t);
fprintf('Time for MEX file: %.3f s\n', t_mex);
fprintf('Speed gain: %.2f\n', t/t_mex);
