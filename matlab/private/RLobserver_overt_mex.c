#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "float.h"

/*
 * RLobserver_overt_mex.c
 *
 *RLOBSERVER_OVERT Responses and log likelihoods for overt RL observer.
 *
 * ================ INPUT VARIABLES ====================
 * PARAMETERS: learning rate and adjustment noise. [1,2] (double)
 * SIGMA: combined category and sensory noise. [scalar] (double)
 * DMU: distance between category means. [scalar] (double)
 * X: matrix of noisy measurements. [Nt,Ns] (double)
 * P_INITIAL: initial probability. [scalar] (double)
 * Z_RESP: subject's criterion responses. [Nt,1] (double)
 * SCORE: Trial feedback. [Nt,1] (double)
 * 
 * ================ OUTPUT VARIABLES ==================
 * MODEL_RESP: model criterion responses, per trial/sample. [Nt,Ns] (double)
 * LOG_PZ: log likelihood, per trial/sample. [Nt,Ns] (double)
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 05-Oct-2016 with MEXXER v0.2 
 * (https://github.com/lacerbi/mexxer).
 */

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

#define M_PI 3.14159265358979323846

void RLobserver_overt( double *model_resp, double *log_Pz, double *parameters, double sigma, double dmu, double *X, double p_initial, double *z_resp, double *score, mwSize Ns, mwSize Nt )
{
    mwSize qq, t;
    double Gamma, m0, a1, a2, *z0;
    
    z0 = z_resp;
    
    Gamma = p_initial/(1-p_initial);
    
    m0 = sigma*sigma * log(Gamma) / dmu;
    a1 = -0.5*log(2*M_PI*parameters[1]);
    a2 = -0.5/(parameters[1]*parameters[1]);
	
	/* Write your main calculations here... */
	for (qq = 0; qq < Ns; qq++) {        
        z_resp = z0;
        
        *(model_resp) = m0;
        *(log_Pz++) = a1 + a2 * (*z_resp-*model_resp)*(*z_resp-*model_resp);
        model_resp++;   z_resp++;   X++;
        
        for (t = 1; t < Nt; t++) {
            if (score[t-1] == 0.) {
                *model_resp = *(model_resp-1) + parameters[0]*( *X - *(model_resp-1) );                
            }
            else {
                *model_resp = *(model_resp-1);
            }
            *(log_Pz++) = a1 + a2 * (*z_resp-*model_resp)*(*z_resp-*model_resp);
            model_resp++;   z_resp++;   X++;
        }
    }

    /*
for qq = 1:Ns
    x = X(:,qq); % vector of noisy measurements
    model_resp(1,qq) = sigma^2 * log(Gamma(1)) / dmu;
    log_Pz(1,qq) = -0.5*log(2*pi*parameters(2)) - 0.5*((z_resp(1)-model_resp(1,qq))./parameters(2)).^2;
    for t = 2:n
        if score(t) == 0
            model_resp(t,qq) = model_resp(t-1,qq) + parameters(1)*(x(t-1)-model_resp(t-1,qq));
        else
            model_resp(t,qq) = model_resp(t-1,qq);
        end
        log_Pz(t,qq) = -0.5*log(2*pi*parameters(2)) - 0.5*((z_resp(t)-model_resp(t,qq))./parameters(2)).^2;
    end
end   
     */ 
    
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *model_resp, *log_Pz, *parameters, sigma, dmu, *X, p_initial, *z_resp, *score;
#if ( ARGSCHECK==0 )
	mwSize *dims_X, *dims_z_resp, *dims_score;
#else /* ( ARGSCHECK!=0 ) */ 
	mwSize *dims_parameters, *dims_X, *dims_z_resp, *dims_score;
#endif /* ( ARGSCHECK!=0 ) */ 
	mwSize Ns, Nt;

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
	   within an if statement, because it will never get to the else
	   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks
	   you out of the MEX-file) */
	if ( nrhs<7 || nrhs>7 )
		mexErrMsgIdAndTxt( "MATLAB:RLobserver_overt:invalidNumInputs",
			"Seven inputs required.");
	if ( nlhs<2 || nlhs>2 )
		mexErrMsgIdAndTxt( "MATLAB:RLobserver_overt:invalidNumOutputs",
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

	/* Get sixth input (Z_RESP, Nt-by-1 double) */
	z_resp = (double*) mxGetPr(prhs[5]);
	dims_z_resp = (mwSize*) mxGetDimensions(prhs[5]);

	/* Get seventh input (SCORE, Nt-by-1 double) */
	score = (double*) mxGetPr(prhs[6]);
	dims_score = (mwSize*) mxGetDimensions(prhs[6]);

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */
#if ( ARGSCHECK==1 )
		dims_parameters = (mwSize*) mxGetDimensions(prhs[0]);
		if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
				mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:parametersNotReal", "Input PARAMETERS must be real.");
		if ( dims_parameters[0] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:parametersWrongSize", "The first dimension of input PARAMETERS has the wrong size (should be 1).");
		if ( dims_parameters[1] != ((mwSize) (2)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:parametersWrongSize", "The second dimension of input PARAMETERS has the wrong size (should be 2).");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetN(prhs[1])*mxGetM(prhs[1])!=1) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:sigmaNotScalar", "Input SIGMA must be a scalar.");

		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || (mxGetN(prhs[2])*mxGetM(prhs[2])!=1) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:dmuNotScalar", "Input DMU must be a scalar.");

		if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )
				mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:XNotReal", "Input X must be real.");
		if ( dims_X[0] != ((mwSize) (Nt)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:XWrongSize", "The first dimension of input X has the wrong size (should be Nt).");
		if ( dims_X[1] != ((mwSize) (Ns)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:XWrongSize", "The second dimension of input X has the wrong size (should be Ns).");

		if ( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || (mxGetN(prhs[4])*mxGetM(prhs[4])!=1) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:p_initialNotScalar", "Input P_INITIAL must be a scalar.");

		if ( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) )
				mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:z_respNotReal", "Input Z_RESP must be real.");
		if ( dims_z_resp[0] != ((mwSize) (Nt)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:z_respWrongSize", "The first dimension of input Z_RESP has the wrong size (should be Nt).");
		if ( dims_z_resp[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:z_respWrongSize", "The second dimension of input Z_RESP has the wrong size (should be 1).");

		if ( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) )
				mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:scoreNotReal", "Input SCORE must be real.");
		if ( dims_score[0] != ((mwSize) (Nt)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:scoreWrongSize", "The first dimension of input SCORE has the wrong size (should be Nt).");
		if ( dims_score[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:RLobserver_overt:scoreWrongSize", "The second dimension of input SCORE has the wrong size (should be 1).");
#endif /* ( ARGSCHECK==1 ) */ 

	/* Pointer to first output (MODEL_RESP, Nt-by-Ns double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (Nt), (mwSize) (Ns), mxREAL);
	model_resp = mxGetPr(plhs[0]);

	/* Pointer to second output (LOG_PZ, Nt-by-Ns double) */
	plhs[1] = mxCreateDoubleMatrix((mwSize) (Nt), (mwSize) (Ns), mxREAL);
	log_Pz = mxGetPr(plhs[1]);

	/* Call the C subroutine */
	RLobserver_overt(model_resp, log_Pz, parameters, sigma, dmu, X, p_initial, z_resp, score, Ns, Nt);

}
