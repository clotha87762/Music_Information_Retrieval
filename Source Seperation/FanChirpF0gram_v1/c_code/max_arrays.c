/* returns the maximum of a series of consecutive concatenated arrays */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray *px;
    mxArray *pArraylengths;
    
    if(nrhs == 2 ){
        /* array values*/
        px = (mxArray*)prhs[0]; 
        /* lengths of the array */
        pArraylengths = (mxArray*)prhs[1]; 
        if(mxIsClass(px, "double")){
            double *yR;
            double *pValXr;
            
            double *pCant_parciales_en_freq;
            int cant_freqs;
            int i;
            int j;
            
            pValXr = mxGetPr(px);
            pCant_parciales_en_freq = mxGetPr(pArraylengths);
            
            cant_freqs = mxGetN(pArraylengths)*mxGetM(pArraylengths);
            
            plhs[0] = mxCreateDoubleMatrix( 1,cant_freqs, mxREAL); 
            yR = mxGetPr(plhs[0]);
            
            for(i = cant_freqs ; i ; i--){
                double max_val = 0;
                int no_parciales = *(pCant_parciales_en_freq++);
                double new_val_min = 0,new_val_max = 0;
                double new_val;
                for(j = no_parciales ; j ;j--){
                  max_val = (max_val< (new_val = *(pValXr++)))?new_val:max_val;
                }
                *(yR++) = max_val;
            }               
        }
        
        if(mxIsClass(px, "single")){
            float *yR;
            float *pValXr;
            double *pCant_parciales_en_freq;
            int cant_freqs;
            int i;
            int j;
            
            pValXr = (float*)mxGetPr(px);
            pCant_parciales_en_freq = mxGetPr(pArraylengths);
            
            cant_freqs = mxGetN(pArraylengths)*mxGetM(pArraylengths);
            plhs[0] = mxCreateNumericMatrix(1, cant_freqs,mxSINGLE_CLASS, mxREAL); /*mxReal is our data-type*/
            yR = (float*)mxGetPr(plhs[0]);
            
            for(i = cant_freqs ; i ; i--){
                double  max_val = 0;
                int no_parciales = *(pCant_parciales_en_freq++);
                
                double new_val;
                for(j = no_parciales ; j ;j--){
                    max_val = (max_val< (new_val = *(pValXr++)))?new_val:max_val;
                }
                *(yR++) = (float)max_val;
            }
        }
    }else{
        printf("max_arrays\nIncorrect input/output parameters.\n\n");
        printf("See help for usage\n");
    }
    return;
}
