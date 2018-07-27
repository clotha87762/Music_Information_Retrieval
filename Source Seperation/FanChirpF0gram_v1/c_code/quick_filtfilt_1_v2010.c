/*You can include any C libraries that you normally use*/
#include "math.h"
#include "stdlib.h"
#include "mex.h"  
#include "stdio.h"

filtfilt_2by1a(double* apIn, double *apOut, long length, double b1, double b2, double a2){
        double *pIn;
        double * pOut;
        double x_n,x_n_1,y_n_1;
        int i;
        x_n_1 = apIn[5];
        y_n_1 = x_n_1;
        for(i=4;i>=1;i--){     
            x_n = apIn[i];/*antisimétrico! 2*pValXr[0]-*/
            apOut[i] = (y_n_1 = b1*x_n_1+b2*x_n - a2*y_n_1);
            x_n_1 = x_n;
        }
        pIn = apIn;
        pOut = apOut;
        for(i=0;i<length-1;i++){     
            x_n = *(pIn++);/*pValXr[i];*/
            *(pOut++) = (y_n_1 = b1*x_n_1+b2*x_n - a2*y_n_1);/*outArrayR[i] = */
            x_n_1 = x_n;
        }
        pOut = apOut+length-1;
        y_n_1 = (x_n_1 = *pOut);
        for(i=length-1;i>=0;i--){     
            x_n = *pOut;/*outArrayR[i];*pOut; */
                   /* printf("%f ",x_n_1);*/
            *(pOut--) = (y_n_1 = b1*x_n_1+b2*x_n - a2*y_n_1);/*outArrayR[i] = */
            x_n_1 = x_n;
        }
}

/* Si b1== b2, se ahorra un producto!!! */
filtfilt_1by1a(double* apIn, double *apOut, long length, double b, double a2){
        double *pIn;
        double * pOut;
        double x_n,x_n_1,y_n_1;
        int i;
       /* x_n_1 = apIn[5];
        y_n_1 = x_n_1;
        for(i=4;i>=1;i--){     
            x_n = apIn[i];
            apOut[i] = (y_n_1 = b*(x_n-x_n_1) - a2*y_n_1);
            x_n_1 = x_n;
        }*/
        x_n_1 = 0;
        y_n_1 = 0;
        
        pIn = apIn;
        pOut = apOut;
   
/*          __asm{
             prefetchnta [apIn + 4096]
          }*/
        for(i=length; i ;i--){     

            *(pOut++) = (y_n_1 = (b*((x_n = *(pIn++))+x_n_1) - a2*y_n_1));/*outArrayR[i] = */
            x_n_1 = x_n;
        }
        pOut = apOut+length-2;

          /*__asm{
             prefetchnta [apIn]
          }*/

        for(i=length-1;i;i--){     
            *(pOut--) = (y_n_1 = b*((x_n = *pOut)+x_n_1) - a2*y_n_1);/*outArrayR[i] = */
            x_n_1 = x_n;
        }
}

/* Si b1== -b2, se ahorra un producto!!! */
filtfilt_1_by1a(double* apIn, double *apOut, long length, double b, double a2){
        double *pIn;
        double * pOut;
        double x_n,x_n_1,y_n_1;
        int i;
       /* x_n_1 = apIn[5];
        y_n_1 = x_n_1;
        for(i=4;i>=1;i--){     
            x_n = apIn[i];
            apOut[i] = (y_n_1 = b*(x_n-x_n_1) - a2*y_n_1);
            x_n_1 = x_n;
        }*/
        x_n_1 = 0;
        y_n_1 = 0;
        
        pIn = apIn;
        pOut = apOut;
   
/*          __asm{
             prefetchnta [apIn + 4096]
          }*/
        for(i=length; i ;i--){     

            *(pOut++) = (y_n_1 = (b*((x_n = *(pIn++))-x_n_1) - a2*y_n_1));/* outArrayR[i] = */
            x_n_1 = x_n;
        }
        pOut = apOut+length-2;

          /*__asm{
             prefetchnta [apIn]
          }*/

        for(i=length-1;i;i--){     
            *(pOut--) = (y_n_1 = b*((x_n = *pOut)-x_n_1) - a2*y_n_1);/*outArrayR[i] = */
            x_n_1 = x_n;
        }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double x_n,x_n_1;
    double y_n,y_n_1;

    /*Declarations*/
    mxArray *pB;
    mxArray *pA;
    mxArray *px;

    double b1,b2,a2;
    double *pValBr;
    double *pValAr;
    double *pValXr;
    double *pValBi;
    double *pValAi;
    double *pValXi;
    int largoB,largoA,largox;
    
    double *outArrayR;    
        double *outArrayI;  
    
    int i;
    
    if(nrhs == 3){
    /*Copy input pointer x*/
        pB = (mxArray *)prhs[0];
        pA = (mxArray *)prhs[1];
        px = (mxArray *)prhs[2];
        
    /*Get matrix x*/
        pValBr = mxGetPr(pB);
        b1 = pValBr[0];
        if(mxGetN(pB)*mxGetM(pB)>1){
            b2 = pValBr[1];
        }else{
            b2 = b1;
        }
 
        pValAr = mxGetPr(pA);
        a2 = pValAr[1];
        pValXr = mxGetPr(px);
        
        largox = mxGetN(px)*mxGetM(px);
        if(1){/*mxIsClass(pB, "double") && mxIsClass(pB, "double") && mxIsClass(pB, "double")){*/
            /*Salida*/
            if(nlhs==1){
                int i;
                double *pOut;
                double *pIn;
                if(mxIsComplex (prhs[2])){
                       plhs[0] = mxCreateDoubleMatrix(mxGetM(px), mxGetN(px), mxCOMPLEX); /*mxReal is our data-type*/                    
                       outArrayI = mxGetPi(plhs[0]);   
                       pValXi = mxGetPi(px);
                       
                }else{
                      plhs[0] = mxCreateDoubleMatrix(mxGetM(px), mxGetN(px), mxREAL); /*mxReal is our data-type*/                    
                          
                }                
                outArrayR = mxGetPr(plhs[0]);
                if(mxGetN(pB)*mxGetM(pB)>1 || b1 != b2 || b1 != -b2){
                    filtfilt_2by1a(pValXr, outArrayR, largox, b1, b2, a2);
                    if(mxIsComplex (prhs[2])){
                        filtfilt_2by1a(pValXi, outArrayI, largox, b1, b2, a2);
                    }    
                }else{
                    if(b1 == b2){
                        filtfilt_1by1a(pValXr, outArrayR, largox, b1, a2);
                        if(mxIsComplex (prhs[2])){
                            filtfilt_1by1a(pValXi, outArrayI, largox, b1, a2);
                        }            
                    }else{
                        filtfilt_1_by1a(pValXr, outArrayR, largox, b1, a2);
                        if(mxIsComplex (prhs[2])){
                            filtfilt_1_by1a(pValXi, outArrayI, largox, b1, a2);
                        }   
                    }
                }                   
               }   
            }
        }else{
        printf("Not enough input arguments...");
    }
    return;
}

