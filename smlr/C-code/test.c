
#if !defined(_WIN32)
#define dgemm dgemm_
#endif

#include <math.h>
#include<stdio.h>
#include "mex.h"
#include "blas.h"
# include<stdlib.h>


void AB (int Ar, int Ac, int Br, int Bc, double* A, double* B, double* C);
void ATB (int Ar, int Ac, int Br, int Bc, double* A, double* B, double* C);
void ABT (int Ar, int Ac, int Br, int Bc, double* A, double* B, double* C);
void sum(mwSize r, mwSize c, double*in, double* out);
void subtract_row (mwSize r1, mwSize c1, double* A, mwSize ir1,
        mwSize r2, mwSize c2, double* B, mwSize ir2, double* C);
void sum2(mwSize r, mwSize c, double*in, double* out);
void power_comp(mwSize r, mwSize c, double*in, double n);

double softmax(double a, double b);
void init(mwSize r, mwSize c, double*in, double scalar);
void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[] )

{
    double *A, *B, *C, *D;

    mwSize k, m, n,m1,n1;
    double* MAXITER;
     /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","2 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","1 output required.");
    }
    
    A = mxGetPr(prhs[0]); /* pointer to first input matrix */
    B = mxGetPr(prhs[1]); /* pointer to second input matrix */
    

    
    /* dimensions of input matrices 
     X [ dxn] Y[mxn]
     A [mxk] B[m,n]*/
    m = mxGetM(prhs[0]);  
    n = mxGetN(prhs[0]);
   m1 = mxGetM(prhs[1]);
   n1 = mxGetN(prhs[1]);
    
   /*I have to check the dimension */
//     plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
//     C = mxGetPr(plhs[0]);
//     subtract_row (m, n, A, 0, m1, n1, B, 1, C);
    
//     plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
//     C = mxGetPr(plhs[0]);
//     memcpy(C,(A+m),m*sizeof(double));
//     C[0]=-1e-10;
   
//    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
//    C = mxGetPr(plhs[0]);
//    memcpy(C,A,m*n*sizeof(double));
//    //sum2(m,n,A,C);   
//    power_comp(m,n,C,2);
   
   
//    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
//    C = mxGetPr(plhs[0]);   
//    C[0]=softmax(A[0],A[1]);

   
      plhs[0] = mxCreateDoubleMatrix(n, n1, mxREAL);
   C = mxGetPr(plhs[0]);  
   D =mxMalloc(n*n1*sizeof(double));
   ATB(m,n,m1,n1,A,B,D);
   memcpy(C,D,n*n1*sizeof(double));


    return;   
}


void sum(mwSize r, mwSize c, double*in, double* out)
{
    mwSize ir, ic;
    
    for (ic=0; ic<c; ic++){
        out[ic]=0;
        for (ir=0; ir<r; ir++){            
            out[ic] += in[ir+ic*r];                
            
        }
    }
    
}


/* This routine performs a dgemm operation
 *  C := C + A * B
*/ 
void AB (int Ar, int Ac, int Br, int Bc, double* A, double* B, double* C)
{
  /* A[MxK] B[KXN]==> AB*/  
  char TRANSA = 'N';
  char TRANSB = 'N';
  double ALPHA = 1.;
  double BETA = 0;  
  dgemm(&TRANSA, &TRANSB, &Ar, &Bc, &Ac, &ALPHA, A, &Ar, B, &Br, &BETA, C, &Ar);
} 
/* This routine performs a dgemm operation
 *  C := C + A' * B
*/ 
void ATB (int Ar, int Ac, int Br, int Bc, double* A, double* B, double* C)
{
    /* A[MxK] B[MXN]==> ATB*/
  char TRANSA = 'T';
  char TRANSB = 'N';
  double ALPHA = 1.;
  double BETA = 0;    
init(Ac, Bc, C, 0.0);
  
  /*DGEMM(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);;*/
  dgemm(&TRANSA, &TRANSB, &Ac, &Bc, &Ar, &ALPHA, A, &Ac, B, &Br, &BETA, C, &Ac);  
} 


/* This routine performs a dgemm operation
 *  C := C + A' * B
*/ 
void ABT (int Ar, int Ac, int Br, int Bc, double* A, double* B, double* C)
{
    /* A[MxN] B[KxN]==> ATB*/
  char TRANSA = 'N';
  char TRANSB = 'T';
  double ALPHA = 1.;
  double BETA = 0;    

  
  /*DGEMM(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);;*/
  dgemm(&TRANSA, &TRANSB, &Ar, &Br, &Ac, &ALPHA, A, &Ar, B, &Bc, &BETA, C, &Ar);  
} 


void subtract_row (mwSize r1, mwSize c1, double* A, mwSize ir1,
        mwSize r2, mwSize c2, double* B, mwSize ir2, double* C){
    
    mwSize ic;
    for (ic=0; ic<c1; ic++){
        C[ic] = A[ir1+r1*ic] - B[ir2+r2*ic];
    }
}
/* sum across column*/
void sum2(mwSize r, mwSize c, double*in, double* out)
{
    mwSize ir, ic;    
    for (ir=0; ir<r; ir++){
        out[ir]=0;
        for (ic=0; ic<c; ic++){
            out[ir] += in[ir+ic*r];                
        }
    }
    
}    
    
    
void power_comp(mwSize r, mwSize c, double*in, double n)
{
    mwSize ir, ic;
    
    for (ic=0; ic<c; ic++){
        for (ir=0; ir<r; ir++){            
            in[ir+ic*r] = pow(in[ir+ic*r],n);                
        }
    }
}    
    
double softmax(double a, double b){
    double out, tmp1, tmp2;

    tmp1 = fabs(a)-b;
    // max([0 abs(a)-b]);
    if (tmp1<0){
        tmp1 =0;
    }
    // sign function
    if(a>0){
        tmp2=1;
    }else{
        if(a<0){
            tmp2=-1;
        }else{
        tmp2=0;
        }
    }
   
    out = tmp1*tmp2;
    return out;
}    
void init(mwSize r, mwSize c, double*in, double scalar)
{
    mwSize ir, ic;
    
    for (ic=0; ic<c; ic++){
        for (ir=0; ir<r; ir++){            
            in[ir+ic*r] = scalar;                
        }
    }
    
} 