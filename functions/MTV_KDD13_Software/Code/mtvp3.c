#include "mtv.h"

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* input =          mxGetPr(prhs[0]);
    const mwSize *dim =      mxGetDimensions(prhs[0]);
    double lam2 =            mxGetScalar(prhs[1]); 
    double rho =             mxGetScalar(prhs[2]);  
    size_t maxIter =         (size_t) mxGetScalar(prhs[3]);    
    double* tol =            mxGetPr(prhs[4]);    
    int display =            (int) mxGetScalar(prhs[5]);
    double *zr, *output;
    size_t numofdim =  (size_t) mxGetNumberOfDimensions(prhs[0]);   
    size_t m, n, d;
    if(numofdim == 3){
       m = (size_t) dim[0];
       n = (size_t) dim[1];
       d = (size_t) dim[2];
       printf("The dimensions are %d x %d x %d\n", (int) m, (int) n, (int) d);
    }else{
        printf("The dimension should be 3...%d\n", (int) numofdim);
        return;
    }
    
    plhs[0] = mxCreateDoubleMatrix(m, n*d, mxREAL);
    mxSetDimensions(plhs[0], dim, 3);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    output = mxGetPr(plhs[0]);
    zr = mxGetPr(plhs[1]);
    memcpy(output, input, m * n * d * sizeof(double));
    zr[0] = (double)mtvp3(output, input, m, n, d, lam2, rho, maxIter,tol[0], tol[1],display);
}















    //double* input =          mxGetPr(prhs[0]);
    //const mwSize *dim =      mxGetDimensions(prhs[0]);
//     double lam2 =            mxGetScalar(prhs[1]); 
//     double rho =             mxGetScalar(prhs[2]);  
//     size_t maxIter =         (size_t) mxGetScalar(prhs[3]);    
//     double* tol =            mxGetPr(prhs[4]);    
//     int display =            (int) mxGetScalar(prhs[5]);
   // double *zr, *output; 
   // size_t numofdim =  1;
            //(size_t) mxGetNumberOfDimensions(prhs[0]);
   // size_t m, n, d;
    //if(numofdim == 3){
       //m = (size_t) dim[0];
       //n = (size_t) dim[1];
       //d = (size_t) dim[2];
   // }else{
    //    printf("The dimension should be 3...\n");
    //}

    //plhs[0] = mxCreateDoubleMatrix(2, 3, mxREAL);
    //mxSetDimensions(plhs[0], dim, 3);
    //plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    //output = mxGetPr(plhs[0]);
    
    //zr = mxGetPr(plhs[1]);
    //memcpy(output, input, m * n * d * sizeof(double));
    //zr[0] = (double)mtvp3(output, input, m, n, d, lam2, rho, maxIter,tol[0], tol[1],display);