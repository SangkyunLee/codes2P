#include <math.h>
#include <stdio.h>
#include "mex.h"


static void eucl_dist( int m, int n, double *input,double *output )
{
    int i,j,k;

	for(i=0;i<n; i++){
		for (j=0; j<=i; j++){
			for (k=0; k<m; k++){
				/*output[i*n+j]+=(input[k*n+i]-input[k*n+j])*(input[k*n+i]-input[k*n+j]);*/
				output[i*n+j]+=(input[k+i*m]-input[k+j*m])*(input[k+i*m]-input[k+j*m]);
			}
		}
	}
    return;
}

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[] )

{
	int m, n;
	int p;
	double *input,*output;
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);



    /* Create an mxArray for the output data */
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);



    /* Retrieve the input data */
    input = mxGetPr(prhs[0]);

    /* Create a pointer to the output data */
    output = mxGetPr(plhs[0]);
	eucl_dist(m,n,input,output);

    return;

}


