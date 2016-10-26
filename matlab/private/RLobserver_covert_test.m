%RLOBSERVER_COVERT_TEST
%   Test script for MEX-file RLobserver_covert.
%
%   Template MATLAB code generated on 25-Oct-2016 with MEXXER v0.2 
%   (https://github.com/lacerbi/mexxer).

TolErr = sqrt(eps);	% Maximum error tolerance per array element

% Define array sizes (choose reasonable values)
Ns = 1000;
Nt = 800;

% Randomly initialize input variables
% (or write here alternative initializations)
parameters = rand([1,2]);	%PARAMETERS: sensory noise and learning rate.
sigma = 10*rand([1,1]);	%SIGMA: combined category and sensory noise.
dmu = 10*rand([1,1]);	%DMU: distance between category means.
X = 10*rand([Nt,Ns]);	%X: matrix of noisy measurements.
p_initial = rand([1,1]);	%P_INITIAL: initial probability.
resp_obs = randi(2,[Nt,1]);	%RESP_OBS: subject's categorization responses.
score = randi(2,[Nt,1])-1;	%SCORE: Trial feedback.

fprintf('==========================\n');
fprintf('Testing RLobserver_covert:\n');
fprintf('==========================\n');

% Call MATLAB and MEX functions
tic; [model_resp,log_P] = RLobserver_covert(parameters,sigma,dmu,X,p_initial,resp_obs,score); t = toc;
tic; [model_resp_mex,log_P_mex] = RLobserver_covert_mex(parameters,sigma,dmu,X,p_initial,resp_obs,score); t_mex = toc;

% Correctness check
model_resp_err = sum(abs(model_resp(:)-model_resp_mex(:)));
fprintf('Total error (model_resp): %g\n', model_resp_err);
if model_resp_err > TolErr*numel(model_resp);
	error('mexxer:tooLargeError','Correctness check failed. Error too large in model_resp.');
end
log_P_err = sum(abs(log_P(:)-log_P_mex(:)));
fprintf('Total error (log_P): %g\n', log_P_err);
if log_P_err > TolErr*numel(log_P);
	error('mexxer:tooLargeError','Correctness check failed. Error too large in log_P.');
end

% Runtime analysis
fprintf('Time for MATLAB code: %.3f s\n', t);
fprintf('Time for MEX file: %.3f s\n', t_mex);
fprintf('Speed gain: %.2f\n', t/t_mex);
