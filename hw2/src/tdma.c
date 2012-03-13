#include "matrix.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /* Lower diagonal with a padded 0 */
    double* a = mxGetPr(prhs[0]);
    
    /* Main Diagonal */
    double* b = mxGetPr(prhs[1]);
    
    /* Upper diagonal */
    double* c = mxGetPr(prhs[2]);
    
    /* right hand side of A * x = b */
    double* v = mxGetPr(prhs[3]);
    
    int i;
    int n = mxGetM(prhs[3]);
    double m;
    
    /* loop over the elements */
    for (i = 1; i < n; ++i) {
        
        /* lower divided by middle */
        m = a[i] / b[i - 1];
        
        /* middle minus ai / bim1 times upper */
        b[i] -= m * c[i - 1];
        
        /* sol - m times previous element in solution */
        v[i] -= m * v[i - 1];
    }
        
    /* 
     * Allocate some memory for the result.
     */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    
    /*
     * Get the output vector in a pointer so that we can manipulate it.
     */
    double* x = mxGetPr(plhs[0]);
    
    /* 
     * The last element of the solution is equal to the last element of the 
     * right hand side divided by the last element of the main diagonal.
     */
    x[n - 1] = v[n - 1] / b[n - 1];
    
    /* backsolve */
    for (i = n - 2; i >= 0; --i) {
        x[i] = (v[i] - c[i] * x[i + 1]) / b[i];
    }    
}
