/* Calculates the log(1+abs(\alpha*x)) of the elements of an array*/

#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*Declarations*/
    mxArray *px;
    
    int largox;
    int i;
    
    if(nrhs == 2 && nlhs == 1){
        px = (mxArray*)prhs[0];
        if(mxIsClass(px, "double")){
            double *pValXr;
            double *pValXi;
            
            double alpha;
            double aux;
            double *y;
            double re,im;
            
            pValXr = mxGetPr(px);
            alpha = *((double*)mxGetPr((mxArray*)prhs[1]));
            
            largox = mxGetN(px)*mxGetM(px);
            /*Output*/
            plhs[0] = mxCreateNumericMatrix(mxGetM(px), mxGetN(px),mxDOUBLE_CLASS, mxREAL); /*mxReal is our data-type*/
            y = mxGetPr(plhs[0]);
            if(mxIsComplex(prhs[0])){
                 /* get the pointer to the imaginary part */
                pValXi = mxGetPi(px);
                for(i = largox; i ; i--){
                    re = *(pValXr++);
                    im = *(pValXi++);
                    *(y++) = log10(1.0 + alpha*sqrt(re*re + im*im));
                }
            }else{
                for(i = largox; i ; i--){
                    *(y++) = log10(1.0 + alpha*fabs(*(pValXr++)));
                }
            }
        }
        
    }else{
        printf("log10_1_alpha_abs_c\nIncorrect input/output parameters.\n\n");
        printf("See help\n");
    }
    return;
}
