/* returns the sum of a series of consecutive concatenated arrays */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray *px;
    mxArray *pArraylengths;
    
    if(nrhs == 2){
        px = (mxArray*)prhs[0];             /* array of values to accumuate */
        pArraylengths = (mxArray*)prhs[1];  /* lengths of each array*/
        if(mxIsClass(px, "double")){ 
            double *yR;
            double *pValXr;
            
            int largox;
            double *pCant_parciales_en_freq;
            int cant_freqs;
            int i;
            int j;
            
            pValXr = mxGetPr(px);
            pCant_parciales_en_freq = mxGetPr(pArraylengths);
            
            cant_freqs = mxGetN(pArraylengths)*mxGetM(pArraylengths);
            
            if(mxIsComplex(px)){ /* Manage Real and complex arrays*/
                double *yI;
                double *pValXi = mxGetPi(px);
                plhs[0] = mxCreateDoubleMatrix( 1,cant_freqs, mxCOMPLEX); 
                yR = mxGetPr(plhs[0]);
                yI = mxGetPi(plhs[0]);
                
                /* iterate and accumulate */
                for(i = cant_freqs ; i ; i--){
                    int no_parciales = *(pCant_parciales_en_freq++);
                    double accumR = 0,accumI = 0;
                    for(j = no_parciales ; j ;j--){
                        accumR += *(pValXr++);
                        accumI += *(pValXi++);
                    }
                    *(yR++) = accumR;
                    *(yI++) = accumI;
                }                
            }else{
                plhs[0] = mxCreateDoubleMatrix( 1,cant_freqs, mxREAL); /*mxReal is our data-type*/
                yR = mxGetPr(plhs[0]);
                
                /* iterate and accumulate */
                for(i = cant_freqs ; i ; i--){
                    double accum = 0;
                    int no_parciales = *(pCant_parciales_en_freq++);
                    for(j = no_parciales ; j ;j--){
                        accum += *(pValXr++);
                    }
                    *(yR++) = accum;
                }
            }
        }
        
        if(mxIsClass(px, "single")){
            float *yR;
            float *pValXr;
            int largox;
            double *pCant_parciales_en_freq;
            int cant_freqs;
            int i;
            int j;
            
            pValXr = (float*)mxGetPr(px);
            pCant_parciales_en_freq = mxGetPr(pArraylengths);
            
            cant_freqs = mxGetN(pArraylengths)*mxGetM(pArraylengths);
            
            if(mxIsComplex(px)){
                float *yI;
                float *pValXi = (float*)mxGetPi(px);
                plhs[0] = mxCreateNumericMatrix(1, cant_freqs,mxSINGLE_CLASS, mxCOMPLEX);; /*mxComplex is our data-type*/
                yR = (float*)mxGetPr(plhs[0]);
                yI = (float*)mxGetPi(plhs[0]);
                
                for(i = cant_freqs ; i ; i--){
                    int no_parciales = *(pCant_parciales_en_freq++);
                    float accumR = 0,accumI = 0;
                    for(j = no_parciales ; j ;j--){
                        accumR += *(pValXr++);
                        accumI += *(pValXi++);
                    }
                    *(yR++) = accumR;
                    *(yI++) = accumI;
                }
            }else{
                plhs[0] = mxCreateNumericMatrix(1, cant_freqs,mxSINGLE_CLASS, mxREAL); /*mxReal is our data-type*/
                yR = (float*)mxGetPr(plhs[0]);
                
                for(i = cant_freqs ; i ; i--){
                    double accum = 0;
                    int no_parciales = *(pCant_parciales_en_freq++);
                    
                    for(j = no_parciales ; j ;j--){
                        accum += *(pValXr++);
                    }
                    *(yR++) = accum;
                }
            }
        }
    }else{
        printf("accum_arrays\nIncorrect input/output parameters.\n\n");
        printf("inputs : x, lengths \n");
        printf("outputs : sums \n\n");
  
        printf(" x = [x_1 x_2 ... x_m]  where x_i are arrays \n");
        printf(" of length l_i, and lengths = [l_1 l2 ... l_m] \n\n");
        printf(" length(x) must be equal to sum(lengths)\n");
        printf(" sums(i) = sum(x_i)\n");
        printf(" x must be a real or complex vector of type: single or double\n");
        printf(" lengths must be a real vector of type: double\n\n");
    }
    return;
}
