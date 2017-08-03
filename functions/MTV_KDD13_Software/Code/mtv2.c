#include "mtv.h"

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	/*set up input arguments */
	double* input =          mxGetPr(prhs[0]);
	size_t m=                 (size_t) mxGetM(prhs[0]);
	size_t n=                 (size_t) mxGetN(prhs[0]);
	double lam2 =             mxGetScalar(prhs[1]);	
	double rho =             mxGetScalar(prhs[2]);	
	size_t maxIter =          (size_t) mxGetScalar(prhs[3]);	
	double* tol =              mxGetPr(prhs[4]);	
	int display =              (int) mxGetScalar(prhs[5]);
    double *zr, *output;
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    output = mxGetPr(plhs[0]);
    zr = mxGetPr(plhs[1]);
    memcpy(output, input, m * n * sizeof(double));
	zr[0] = (double)mtv(output, input, m, n, lam2, rho, maxIter,tol[0],tol[1],display);
}