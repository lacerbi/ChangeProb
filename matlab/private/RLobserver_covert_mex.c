#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "float.h"

/*
 * RLobserver_covert_mex.c
 *
 *RLOBSERVER_COVERT Responses and log likelihoods for covert RL observer.
 *
 * ================ INPUT VARIABLES ====================
 * PARAMETERS: sensory noise and learning rate. [1,2] (double)
 * SIGMA: combined category and sensory noise. [scalar] (double)
 * DMU: distance between category means. [scalar] (double)
 * X: matrix of noisy measurements. [Nt,Ns] (double)
 * P_INITIAL: initial probability. [scalar] (double)
 * RESP_OBS: subject's categorization responses. [Nt,1] (double)
 * SCORE: Trial feedback. [Nt,1] (double)
 * 
 * ================ OUTPUT VARIABLES ==================
 * MODEL_RESP: model categorization responses probability, per trial/sample. [Nt,Ns] (double)
 * LOG_P: log likelihood, per trial/sample. [Nt,Ns] (double)
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 25-Oct-2016 with MEXXER v0.2 
 * (https://github.com/lacerbi/mexxer).
 */

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

#define M_PI 3.14159265358979323846

void RLobserver_covert( double *model_resp, double *log_P, double *parameters, double sigma, double dmu, double *X, double p_initial, double *resp_obs, double *score, mwSize Ns, mwSize Nt )
{
	
    mwSize qq, t;
    double Gamma, m0, a1, a2, *C0, z_model, lambda;
    double mr0,mr1,logmr1a,logmr1b,logmr0a,logmr0b;
    
    lambda = 0.01; /* Minimum lapse to avoid numerical trouble */
    C0 = resp_obs;
    
    Gamma = p_initial/(1-p_initial);    
    m0 = sigma*sigma * log(Gamma) / dmu;    
	
    /* Precomputed variables to speed up computations */
    mr1 = (1-lambda)*1. + 0.5*lambda;
    mr0 = (1-lambda)*0. + 0.5*lambda;    
    logmr1a = log(mr1);
    logmr1b = log(1.0-mr1);
    logmr0a = log(mr0);
    logmr0b = log(1.0-mr0);
    
	/* Write your main calculations here... */
	for (qq = 0; qq < Ns; qq++) {        
        resp_obs = C0;
        z_model = m0;
        
        if (*(X++) <= z_model) {
            *(model_resp++) = mr1;
            *(log_P++) = (*(resp_obs++) == 1) ? logmr1a : logmr1b;
        }
        else {
            *(model_resp++) = mr0;
            *(log_P++) = (*(resp_obs++) == 1) ? logmr0a : logmr0b;
        }
        
        for (t = 1; t < Nt; t++) {
            if (score[t-1] == 0.) 
                z_model += parameters[1]*( *X - z_model );                
            
            if (*(X++) <= z_model) {
                *(model_resp++) = mr1;
                *(log_P++) = (*(resp_obs++) == 1) ? logmr1a : logmr1b;
            }
            else {
                *(model_resp++) = mr0;
                *(log_P++) = (*(resp_obs++) == 1) ? logmr0a : logmr0b;
            }
            
        }
    }

    /*
for qq = 1:Ns
    x = X(:,qq); % vector of noisy measurements
    Gamma(1) = p_initial/(1-p_initial);
    z_model(1,qq) = sigma^2 * log(Gamma(1)) / dmu;
    if x(1) <= z_model(1,qq)
        model_resp(1,qq) = 1;
    else
        model_resp(1,qq) = 0;
    end
    model_resp(1,qq) = lambda/2 + (1-lambda)*model_resp(1,qq);
    log_P(1,qq) = log(model_resp(1,qq)).*(Chat(1)==1) + log(1-model_resp(1,qq)).*(Chat(1)~=1);
    for t = 2:Nt
        if score(t-1) == 0
            z_model(t,qq) = z_model(t-1,qq) + parameters(2)*(x(t)-z_model(t-1,qq));
        else
            z_model(t,qq) = z_model(t-1,qq);
        end
        if x(t) <= z_model(t,qq)
            model_resp(t,qq) = 1;
        else
            model_resp(t,qq) = 0;
        end
        model_resp(t,qq) = lambda/2 + (1-lambda)*model_resp(t,qq);
        log_P(t,qq) = log(model_resp(t,qq)).*(Chat(t)==1) + log(1-model_resp(t,qq)).*(Chat(t)~=1);
    end
end
     
     */ 
	
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *model_resp, *log_P, *parameters, sigma, dmu, *X, p_initial, *resp_obs, *score;
#if ( ARGSCHECK==0 )
	mwSize *dims_X, *dims_resp_obs, *dims_score;
#else /* ( ARGSCHECK!=0 ) */ 
	mwSize *dims_parameters, *dims_X, *dims_resp_obs, *dims_score;
#endif /* ( ARGSCHECK!=0 ) */ 
	mwSize Ns, Nt;

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
	   within an if statement, because it will never get to the else
	   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks
	   you out of the MEX-file) */
	if ( nrhs<7 || nrhs>7 )
		mexErrMsgIdAndTxt( "MATLAB:RLobserver_covert:invalidNumInputs",
			"Seven inputs required.");
	if ( nlhs<2 || nlhs>2 )
		mexErrMsgIdAndTxt( "MATLAB:RLobserver_covert:invalidNumOutputs",
			"Two outputs required.");

	/* Get first input (PARAMETERS, 1-by-2 double) */
	parameters = (double*) mxGetPr(prhs[0]);

	/* Get second input (SIGMA, scalar double) */
	sigma = (double) mxGetScalar(prhs[1]);

	/* Get third input (DMU, scalar double) */
	dmu = (double) mxGetScalar(prhs[2]);

	/* Get fourth input (X, Nt-by-Ns double) */
	X = (double*) mxGetPr(prhs[3]);
	dims_X = (mwSize*) mxGetDimensions(prhs[3]);
	Nt = dims_X[0];
	Ns = dims_X[1];

	/* Get fifth input (P_INITIAL, scalar double) */
	p_initial = (double) mxGetScalar(prhs[4]);

	/* Get sixth input (RESP_OBS, Nt-by-1 double) */
	resp_obs = (double*) mxGetPr(prhs[5]);
	dims_resp_obs = (mwSize*) mxGetDimensions(prhs[5]);

	/* Get seventh input (SCORE, Nt-by-1 double) */
	score = (double*) mxGetPr(prhs[6]);
	dims_score = (mwSize*) mxGetDimensions(prhs[6]);

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */
#if ( ARGSCHECK==1 )
		dims_parameters = (mwSize*) mxGetDimensions(prhs[0]);
		if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
				mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:parametersNotReal", "Input PARAMETERS must be real.");
		if ( dims_parameters[0] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:parametersWrongSize", "The first dimension of input PARAMETERS has the wrong size (should be 1).");
		if ( dims_parameters[1] != ((mwSize) (2)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:parametersWrongSize", "The second dimension of input PARAMETERS has the wrong size (should be 2).");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetN(prhs[1])*mxGetM(prhs[1])!=1) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:sigmaNotScalar", "Input SIGMA must be a scalar.");

		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || (mxGetN(prhs[2])*mxGetM(prhs[2])!=1) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:dmuNotScalar", "Input DMU must be a scalar.");

		if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )
				mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:XNotReal", "Input X must be real.");
		if ( dims_X[0] != ((mwSize) (Nt)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:XWrongSize", "The first dimension of input X has the wrong size (should be Nt).");
		if ( dims_X[1] != ((mwSize) (Ns)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:XWrongSize", "The second dimension of input X has the wrong size (should be Ns).");

		if ( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || (mxGetN(prhs[4])*mxGetM(prhs[4])!=1) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:p_initialNotScalar", "Input P_INITIAL must be a scalar.");

		if ( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) )
				mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:resp_obsNotReal", "Input RESP_OBS must be real.");
		if ( dims_resp_obs[0] != ((mwSize) (Nt)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:resp_obsWrongSize", "The first dimension of input RESP_OBS has the wrong size (should be Nt).");
		if ( dims_resp_obs[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:resp_obsWrongSize", "The second dimension of input RESP_OBS has the wrong size (should be 1).");

		if ( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) )
				mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:scoreNotReal", "Input SCORE must be real.");
		if ( dims_score[0] != ((mwSize) (Nt)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:scoreWrongSize", "The first dimension of input SCORE has the wrong size (should be Nt).");
		if ( dims_score[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_covert:scoreWrongSize", "The second dimension of input SCORE has the wrong size (should be 1).");
#endif /* ( ARGSCHECK==1 ) */ 

	/* Pointer to first output (MODEL_RESP, Nt-by-Ns double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (Nt), (mwSize) (Ns), mxREAL);
	model_resp = mxGetPr(plhs[0]);

	/* Pointer to second output (LOG_P, Nt-by-Ns double) */
	plhs[1] = mxCreateDoubleMatrix((mwSize) (Nt), (mwSize) (Ns), mxREAL);
	log_P = mxGetPr(plhs[1]);

	/* Call the C subroutine */
	RLobserver_covert(model_resp, log_P, parameters, sigma, dmu, X, p_initial, resp_obs, score, Ns, Nt);

}
