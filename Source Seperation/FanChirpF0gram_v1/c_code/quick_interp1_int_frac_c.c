/* This function performs a linear interpolation a vector 
 * given the integer and fractional part of the
 * interpolating position
 *             y(i,j) = x(pos_int(i,j)+pos_frac(i,j))
 */

/* inputs : x, pos_int (NxM) ,pos_frac (NxM)
 * outputs : y (NxM)*/

#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*Declarations*/
    mxArray *px;
    mxArray *pPosInt;
    mxArray *pPosFrac;
    
    int largox;
    int i;
    
    if(nrhs == 3 && nlhs ==1){
        px = (mxArray*)prhs[0];
        pPosInt = (mxArray*)prhs[1];
        pPosFrac = (mxArray*)prhs[2];
        
        if(mxIsClass(px, "single")){
            float *pValX;
            int *pValPosInt;
            float *pValPosFrac;
            
            float *y;
            
            pValX = (float*)mxGetPr(px);
            pValPosInt = (int*)mxGetPr(pPosInt);
            pValPosFrac = (float*)mxGetPr(pPosFrac);
            
            
            largox = mxGetN(pPosInt)*mxGetM(pPosInt);
            /*Salida*/
            plhs[0] = mxCreateNumericMatrix(mxGetM(pPosInt), mxGetN(pPosInt),mxSINGLE_CLASS, mxREAL); /*mxReal is our data-type*/
            y = (float*)mxGetPr(plhs[0]);
            for(i = 0;i<largox; i++){
                *(y++) =  (pValX[(*pValPosInt)+1]*(*pValPosFrac) + (((float)1.0)-(*(pValPosFrac++)))*pValX[(*(pValPosInt++))]);
            }
        }
        if(mxIsClass(px, "double")){
            double *pValX;
            int *pValPosInt;
            double *pValPosFrac;
            double *y;
            
            pValX = mxGetPr(px);
            pValPosInt = (int*)mxGetPr(pPosInt);
            pValPosFrac = mxGetPr(pPosFrac);
            
            largox = mxGetN(pPosInt)*mxGetM(pPosInt);
            /*Salida*/
            plhs[0] = mxCreateNumericMatrix(mxGetM(pPosInt), mxGetN(pPosInt),mxDOUBLE_CLASS, mxREAL); /*mxReal is our data-type*/
            y = mxGetPr(plhs[0]);
            
            /* pipeline parallelizing*/
            for(i = largox/2; i ; i--){
                double frac,frac_1;
                frac_1 = 1.0 - (frac = *(pValPosFrac++));
                *(y++) =  (pValX[(*pValPosInt)+1]*frac + frac_1*pValX[(*(pValPosInt++))]);
                
                frac_1 = 1.0 - (frac = *(pValPosFrac++));
                *(y++) =  (pValX[(*pValPosInt)+1]*frac + frac_1*pValX[(*(pValPosInt++))]);
            }            
            if(largox & 0x00000001){
                *(y++) =  (pValX[(*pValPosInt)+1]*(*pValPosFrac) + ((1.0)-(*(pValPosFrac++)))*pValX[(*(pValPosInt++))]);
            }
        }
        
    }else{
        printf("quick_interp1_int_frac_c\nIncorrect input/output parameters.\n\n");
        printf("inputs : x, pos_int (NxM) ,pos_frac (NxM) \n");
        printf("outputs : y (NxM) \n\n");
  
        printf("y(i,j) = x(1+pos_int(i,j)+pos_frac(i,j))  (linear interpolation)\n");
        printf("x must be a real vector of type: single or double\n");
        printf("pos_int must be a real matrix of size NxM and type: uint32\n");
        printf("pos_frac must be a real matriz of size NxM and type : single or double\n");
        printf("y is a real matrix of size NxM and type : single or double\n \n");  
    }
    return;
}
