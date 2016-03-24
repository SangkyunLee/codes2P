
#if !defined(_WIN32)
#define dgemm dgemm_
#endif

#include <math.h>
#include<stdio.h>
#include "mex.h"
#include "blas.h"


void AB (int Ar, int Ac, int Br, int Bc, double* A, double* B, double* C);
void ATB (int Ar, int Ac, int Br, int Bc, double* A, double* B, double* C);
void ABT (int Ar, int Ac, int Br, int Bc, double* A, double* B, double* C);
void sum(int r, int c, double* in, double* out);
void sum2(int r, int c, double* in, double* out);
void exp_comp(int r, int c, double *in);
void power_comp(int r, int c, double*in, double n);
void div_dot(int r, int c, double*A, double*B);

        
void subtract_row(int r1, int c1, double* A, int ir1,
        int r2, int c2, double* B, int ir2, double* C);
double softmax(double a, double b);


void addM(int r, int c, double* A, double *B, double scalar);
void addrowvector(int c, double* A, int r2, int c2, double *B, int ir,  double scalar);


double get_LL(int m, int n, double* Y, double* WtX, double* exp_WtX);
double get_Lp(int d, int m, double *W, double LL, double lambda1, double lambda2);
void copyresult(int iter, double* out, double* LL, double* Lp, double* Wdiff);
// 

void init(int r, int c, double *in, double scalar);
void get_Bkkc(int d, int n, double* X, double* Bkkc, int m);
void get_delta_kkc(int d, double*Bkkc, double* delta_kkc,
        double lambda1, double lambda2);


void smlr(double* W0, double* X, double* Y, int d, int m, int n,
        double* W, double* out,double lambda1, double lambda2, int MAXITER);     




void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[] )

{
    double *W0, *X, *Y;
    double *W, *out;
    double *lambda1, *lambda2;
    int d, m, n;
    double* MAXITER;
     /* check for proper number of arguments */
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Six inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two output required.");
    }
    
    W0 = mxGetPr(prhs[0]); /* pointer to first input matrix */
    X = mxGetPr(prhs[1]); /* pointer to second input matrix */
    Y = mxGetPr(prhs[2]); /* pointer to third input matrix */      
    lambda1 = mxGetPr(prhs[3]);
    lambda2 = mxGetPr(prhs[4]);
    MAXITER = mxGetPr(prhs[5]);
    
    /* dimensions of input matrices 
    W[dxm] X[dxn] Y[mxn]  */
    d = mxGetM(prhs[0]);  
    m = mxGetN(prhs[0]);
    n= mxGetN(prhs[1]);
    
   /*I have to check the dimension */
    plhs[0] = mxCreateDoubleMatrix(d, m, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(3, MAXITER[0], mxREAL);
    W = mxGetPr(plhs[0]);
    out = mxGetPr(plhs[1]);
    /*smlr(W0,X,Y,m,n,d,W,out,lambda1[0],lambda2[0],(int)MAXITER[0]);*/
	mexPrintf("%f, %f, %f, %d,%d,%d\n",W0[0],X[0],Y[0],d,m,n);
    smlr(W0,X,Y,d,m,n,W,out,lambda1[0],lambda2[0], (int) MAXITER[0]);
    return;   
}

void smlr(double* W0, double* X, double* Y, int d, int m, int n,
        double* W, double* out,double lambda1, double lambda2, int MAXITER)        
{
    /* dimensions of input matrices 
     X [ dxn] Y[mxn]  W[dxm]*/
    double *Bkkc, *delta_kkc;
    double *WtX, *exp_WtX, *denom;
    double *LL, *Lp, *Wdiff;
    int iter, inx_i, inx_j;
    double *Wi, *WiX, *exp_WiX, *exp_WiX_new, *Prob_i;
    double *Ai, *grad_Wi;
    
    double Wji, Wji_old;
    
    double diff,diff1;
    
    Bkkc = mxMalloc(d*1*sizeof(double)); 
    delta_kkc = mxMalloc(d*1*sizeof(double)); 
    
    //X^2, W^2,     
   // Wp2 = mxMalloc(d*m*sizeof(double));
    
    WtX = mxMalloc(m*n*sizeof(double));
    exp_WtX = mxMalloc(m*n*sizeof(double));
    denom = mxMalloc(1*n*sizeof(double));
    
    
    LL = mxMalloc(1*MAXITER*sizeof(double));
    Lp = mxMalloc(1*MAXITER*sizeof(double));
    Wdiff = mxMalloc(1*MAXITER*sizeof(double));
    
    Wi = mxMalloc(d*1*sizeof(double));
    WiX = mxMalloc(1*n*sizeof(double));
    exp_WiX = mxMalloc(1*n*sizeof(double));
    exp_WiX_new = mxMalloc(1*n*sizeof(double));
    Prob_i = mxMalloc(1*n*sizeof(double));
    
    
    Ai = mxMalloc(1*n*sizeof(double));
    grad_Wi = mxMalloc(d*1*sizeof(double));
   


    /*// constant over iterations
    Bkkc = -1/2*((m-1)/m)*(sum(X.^2,2));     % [d x 1], Bkk=[Bkkc; Bkkc; ... Bkkc];[d(m-1) x 1]
    delta_kkc = -1*lambda2./(Bkkc-lambda1);   % [d x 1] 
     */
    get_Bkkc(d,n, X, Bkkc, m);
    get_delta_kkc(d, Bkkc, delta_kkc, lambda1, lambda2);
    
    
    /*W=W0*/
    memcpy(W, W0, d*m*sizeof(double)); 
    
    /*WtX=W'*X;  %[d x m-1]
    denom = sum(exp(WtX));*/
    ATB(d,m,d,n,W,X,WtX);    
    memcpy(exp_WtX, WtX, m*n*sizeof(double));
    exp_comp(m,n,exp_WtX);    
    sum(m, n, exp_WtX, denom);
    
    /*LL=-Inf*ones(1,MAX_ITER);
    Lp=-Inf*ones(1,MAX_ITER);
    Wdiff=-Inf*ones(1,MAX_ITER);
     */
    init(1, MAXITER, LL, -1e-10);
    memcpy(Lp, LL, MAXITER*sizeof(double));
    memcpy(Wdiff, LL, MAXITER*sizeof(double));
    
    for (iter = 0; iter<MAXITER; iter++){
        Wdiff[iter]=0;
        
        for (inx_i=0; inx_i<m-1; inx_i++){
            
            //Wi
            memcpy(Wi,W+d*inx_i,d*sizeof(double));
            
            
            for (inx_j=0; inx_j<d; inx_j++){
                if (inx_j == 0){
                    /*
                     WiX = W(:,inx_i)'*X;
                     exp_WiX = exp(WiX); 
                     Prob_i = exp_WiX./denom;
                     */
					 
                     ATB(d,1,d,n,Wi,X,WiX);
                     memcpy(exp_WiX, WiX, 1*n*sizeof(double));
                     exp_comp(1,n,exp_WiX);
                     memcpy(Prob_i, exp_WiX, 1*n*sizeof(double));
                     div_dot(1, n, Prob_i, denom);                     
                }
                
                /*Ai = Y(inx_i,:)-Prob_i; Ai[1xn]*/
                subtract_row (m,n,Y, inx_i,1, n, Prob_i, 0, Ai);
                /* grad_Wi = X*Ai'; grad_Wi [ d x1]*/
				
                ABT(d, n,  1, n,X, Ai, grad_Wi);
                
                //Wji_new = (Bkkc(inx_j)*W(inx_j,inx_i) - grad_Wi(inx_j))/(Bkkc(inx_j)-lambda1);
                Wji = (Bkkc[inx_j]*W[inx_j+d*inx_i] -grad_Wi[inx_j])/(Bkkc[inx_j]-lambda1);
                
                //Wji_new = softmax(Wji_new,delta_kkc(inx_j));
                Wji = softmax(Wji, delta_kkc[inx_j]);
                

                Wji_old = W[inx_j+d*inx_i];                
                Wdiff[iter] += fabs(Wji-Wji_old)/((m-1)*d);
                
                //WiX_new = WiX + (Wji_new - Wji_old)*X(inx_j,:);
                /* I overload WiX without defining WiX_new*/
                addrowvector(n, WiX, d, n, X, inx_j,  Wji-Wji_old);
                
                // exp_WiX_new = exp(WiX_new);
                memcpy(exp_WiX_new, WiX, 1*n*sizeof(double));
                exp_comp(1,n,exp_WiX_new);
                
                //denom = denom - exp_WiX + exp_WiX_new;
                addrowvector(n, denom, 1, n, exp_WiX,0,  -1);
                addrowvector(n, denom, 1, n, exp_WiX_new,0,  1);
                
                //update
                /*
                W(inx_j,inx_i)=Wji_new;
                Prob_i = exp_WiX_new./denom;            
                WiX = WiX_new;
                exp_WiX = exp_WiX_new;
                */
                
                W[inx_j+d*inx_i]=Wji;                
                //ATB(d,1,d,n,Wi,X,WiX);                    
                memcpy(Prob_i, exp_WiX_new, 1*n*sizeof(double));
                div_dot(1, n, Prob_i, denom);
                memcpy(exp_WiX, exp_WiX_new, 1*n*sizeof(double));


            }
        }
        // W'X
        ATB(d,m,d,n,W,X,WtX);
        //exp(W'X)
		memcpy(exp_WtX, WtX, m*n*sizeof(double));
        exp_comp(m,n,exp_WtX);
        LL[iter] = get_LL(m, n,Y,WtX, exp_WtX);        
        Lp[iter] = get_Lp(d, m, W, LL[iter], lambda1, lambda2);
        if (iter>0){
            diff = fabs(Lp[iter]-Lp[iter-1]);
            if (iter==1) diff1 = diff;
            if (diff<1e-3 | (diff/diff1)<1e-5){
                mexPrintf("Converged in iteration: %d",iter);
                break;
            }
        }
        
    }
    
    
    copyresult(iter, out, LL, Lp, Wdiff);
//     double *Bkkc, *delta_kkc;
//     double *WtX, *exp_WtX, *denom;
//     double *LL, *Lp, *Wdiff;
//     int iter, inx_i, inx_j;
//     double *Wi, *WiX, *exp_WiX, *exp_WiX_new, *Prob_i;
//     double *Ai, *grad_Wi;
    
    // memory free
    mxFree(Bkkc);
    mxFree(delta_kkc);
    mxFree(WtX);
    mxFree(exp_WtX); 
    mxFree(denom);
    mxFree(LL);
    mxFree(Lp);
    mxFree(Wdiff);    
    mxFree(WiX);
    mxFree(Wi);
    mxFree(exp_WiX);
    mxFree(exp_WiX_new);
    mxFree(Prob_i);
    mxFree(Ai);
    mxFree(grad_Wi);
}



/* This routine performs a dgemm operation
 *  C =  A' * B
*/ 
void ATB (int Ar, int Ac, int Br, int Bc, double* A, double* B, double* C)
{
    /* A[MxK] B[MXN]==> ATB*/
//   char TRANSA = 'T';
//   char TRANSB = 'N';
//   double ALPHA = 1.0;
//   double BETA = 0.0;   
// 
// 
//   /*DGEMM(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);;*/
//   dgemm(&TRANSA, &TRANSB, &Ac, &Bc, &Ar, &ALPHA, A, &Ac, B, &Br, &BETA, C, &Ac);  
    
    int ic1, ic2, ir;
    double tmp;
    for(ic1=0; ic1<Ac; ic1++){
        for(ic2=0; ic2<Bc; ic2++){
            tmp=0;
            for(ir=0; ir<Ar; ir++)
                tmp  += A[ir + Ar*ic1]*B[ir+ic2*Br]; 
            C[ic1+ic2*Ac]=tmp;
        }
    }
//     mexPrintf("\n\n %f \n\n",C[0]);
} 


/* This routine performs a dgemm operation
 *  C =  A * B'
*/ 
void ABT (int Ar, int Ac, int Br, int Bc, double* A, double* B, double* C)
{
    /* A[MxN] B[KxN]==> ATB*/
//   char TRANSA = 'N';
//   char TRANSB = 'T';
//   double ALPHA = 1.0;
//   double BETA = 0.0;   
// 
//   /*DGEMM(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);;*/
//   dgemm(&TRANSA, &TRANSB, &Ar, &Br, &Ac, &ALPHA, A, &Ar, B, &Bc, &BETA, C, &Ar);  
    
    
    int ir1, ir2, ic;
    double tmp;
    for(ir1=0; ir1<Ar; ir1++){
        for(ir2=0; ir2<Br; ir2++){
            tmp=0;
            for(ic=0; ic<Ac; ic++)
                tmp  += A[ir1 + Ar*ic]*B[ir2+ic*Br]; 
            C[ir1+ir2*Ar]=tmp;
        }
    }
} 
  
void sum(int r, int c, double*in, double* out)
{
    int ir, ic;
    
    for (ic=0; ic<c; ic++){
		out[ic]=0;
        for (ir=0; ir<r; ir++){            
            out[ic] += in[ir+ic*r];                
        }
    }
    
}
/* sum across column*/
void sum2(int r, int c, double*in, double* out)
{
    int ir, ic;    
    for (ir=0; ir<r; ir++){
        out[ir]=0;
        for (ic=0; ic<c; ic++){
            out[ir] += in[ir+ic*r];                
        }
    }
    
}    

void exp_comp(int r, int c, double*in)
{
    int ir, ic;
    
    for (ic=0; ic<c; ic++){
        for (ir=0; ir<r; ir++){            
            in[ir+ic*r] = exp(in[ir+ic*r]);                
        }
    }
    
}
void init(int r, int c, double*in, double scalar)
{
    int ir, ic;
    
    for (ic=0; ic<c; ic++){
        for (ir=0; ir<r; ir++){            
            in[ir+ic*r] = scalar;                
        }
    }
    
} 

void div_dot(int r, int c, double*A, double*B)
{
    int ir, ic;
    
    for (ic=0; ic<c; ic++){
        for (ir=0; ir<r; ir++){            
            A[ir+ic*r] = A[ir+ic*r]/B[ir+ic*r];                
        }
    }
    
}

void subtract_row(int r1, int c1, double* A, int ir1,
        int r2, int c2, double* B, int ir2, double* C)
{
    
    int ic;
    for (ic=0; ic<c1; ic++){
        C[ic] = A[ir1+r1*ic] - B[ir2+r2*ic];
    }
}

/* component-wise power*/    
void power_comp(int r, int c, double*in, double n)
{
    int ir, ic;
    
    for (ic=0; ic<c; ic++){
        for (ir=0; ir<r; ir++){            
            in[ir+ic*r] = pow(in[ir+ic*r],n);                
        }
    }
}  

void get_Bkkc(int d,int n, double*X, double* Bkkc, int m)
{
    double *Xp2;
    int ir;
    
    Xp2 = mxMalloc(d*n*sizeof(double)); 
    memcpy(Xp2,X,d*n*sizeof(double)); 
    power_comp(d, n, Xp2, 2);
    sum2(d, n, Xp2, Bkkc);
    
    for (ir=0; ir<d; ir++){
        Bkkc[ir] = -0.5*Bkkc[ir]*((double)m-1.)/((double)m);
    }
    mxFree(Xp2);
}
    
void get_delta_kkc(int d, double*Bkkc, double* delta_kkc, double lambda1, double lambda2)
{
    int ir;    
    for (ir=0; ir<d; ir++){
        delta_kkc[ir] = -1*lambda2/(Bkkc[ir]-lambda1);
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

/* A = A+ scalar*B */
void addM(int r, int c, double* A, double *B,  double scalar){
    int ir, ic;
    
    for (ic=0; ic<c; ic++){
        for (ir=0; ir<r; ir++){            
            A[ir+ic*r] += B[ir+ic*r]*scalar;                
        }
    }
}
/* a = a+ scalar*B(ir,:) */
void addrowvector(int c, double* A, int r2, int c2, double *B, int ir,  double scalar){
    int ic;    
    for (ic=0; ic<c; ic++){            
        A[ic] += B[ir+ic*r2]*scalar;                
    }
    
}

double get_LL(int m, int n, double* Y, double* WtX, double* exp_WtX)
{
     /* dimensions of input matrices 
     X [ dxn] Y[mxn]  W[dxm]*/
    double* K;
    double LL;
    int ir, ic;
    
    K = mxMalloc(1*n*sizeof(double)); 
    sum(m,n, exp_WtX, K);
    LL = 0;    
    for (ic=0; ic<n; ic++){
        for (ir=0; ir<m; ir++){ 
            LL += WtX[ir+ic*m]*Y[ir+ic*m];            
        }        
    }
    for (ic=0; ic<n; ic++){       
        LL -= log(K[ic]);
    }
    return LL;
}

double get_Lp(int d, int m, double *W, double LL, double lambda1, double lambda2)
{
    double Lp;
    int ir, ic;
    Lp=LL;
    for (ic=0; ic<m; ic++){
        for (ir=0; ir<d; ir++){ 
            Lp -= lambda1*pow(W[ir+ic*d],2);
            Lp -= lambda2*fabs(W[ir+ic*d]);
        }        
    }
    return Lp;
}

void copyresult(int iter, double* out, double* LL, double* Lp, double* Wdiff){
    
    int i;    
    for (i=0; i<=iter; i++){
        out[0 + i*3]=LL[i];
        out[1 + i*3]=Lp[i];
        out[2 + i*3]=Wdiff[i];
    }
}
        
    
    
    

