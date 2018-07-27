/*Applies the iir_cqt filtering to a stft vector given the designed poles*/

#include "mex.h"
#include "stdio.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*Declarations*/
    mxArray *px;
    mxArray *pPoles;
    long prepad;
    long pospad;
    
    
    int largox;
    int i;
    
    if(nrhs == 4){
        
        px = (mxArray*)prhs[0];
        pPoles = (mxArray*)prhs[1];
        prepad = (long)(*((double*)mxGetPr((mxArray*)prhs[2])));
        pospad = (long)(*((double*)mxGetPr((mxArray*)prhs[3])));
        if(mxIsClass(px, "double")){/* && mxIsClass(pB, "double") && mxIsClass(pB, "double")){*/
            double *pValXr;
            double *pValXi;
            double *pValPoles;
            double *yr;
            double *yi;
            double *yout_r;
            double *yout_i;
            mxArray * y_aux;
            long largo_y;
            double aux_pole;
            if(mxIsComplex(px)){
                pValXr = (double*)mxGetPr(px);
                pValXi = (double*)mxGetPi(px);
                pValPoles = (double*)mxGetPr(pPoles);
                
                largox = mxGetN(pPoles)*mxGetM(pPoles);
                largo_y = largox-prepad-pospad;
                
            /*Salida*/
            /* plhs[0] = mxCreateNumericMatrix(mxGetM(pPoles), mxGetN(pPoles),mxDOUBLE_CLASS, mxREAL); */
                y_aux = mxCreateNumericMatrix(mxGetM(pPoles), mxGetN(pPoles),mxDOUBLE_CLASS, mxCOMPLEX);
                plhs[0] = mxCreateNumericMatrix(1,largo_y,mxDOUBLE_CLASS, mxCOMPLEX);
                
                yr = (double*)mxGetPr(y_aux);
                yi = (double*)mxGetPi(y_aux);
                
                yr[prepad] = pValXr[prepad];
                yi[prepad] = -pValXi[prepad];
                
             /*Pre-pad filtering*/
                for(i = prepad-1; i>=0 ; i--){

                    yr[i] = (pValXr[i+1] + pValXr[i]) + yr[i+1]*( aux_pole = pValPoles[prepad-1-i]) ;
                    yi[i] = -(pValXi[i+1] + pValXi[i]) + yi[i+1]* aux_pole;
                }
             /* Move the pointer to have indexes aligned to output indexes */
                pValPoles = pValPoles + prepad;
                
             /* forward filtering */
                for(i = 1; i<largox-prepad ; i++){
                    yr[i] = (pValXr[i-1] + pValXr[i]) + yr[i-1]* ( aux_pole = pValPoles[i]);
                    yr[i-1] = yr[i-1] + yr[i];
                    yi[i] = (pValXi[i-1] + pValXi[i]) + yi[i-1]* aux_pole;
                    yi[i-1] += yi[i];
                }
                
             /* backward filtering */
                for(i = largox-prepad-2; i>=0 ; i--){
                    yr[i] += yr[i+1]*( aux_pole = pValPoles[i]);
                    yi[i] += yi[i+1]* aux_pole;
                }
                yout_r = (double*)mxGetPr(plhs[0]);
                yout_i = (double*)mxGetPi(plhs[0]);
                
            /* Gain correction  & output storage */
                for(i = 0; i<largo_y ; i++){
                    aux_pole = (1.0-(*(pValPoles++)))/2.0;
                    *(yout_r++) = (*(yr++))*aux_pole*aux_pole;/* /((double)1024)/2.0; */
                    *(yout_i++) = (*(yi++))*aux_pole*aux_pole;/* /((double)1024)/2.0; */
                }
            }else{
                pValXr = (double*)mxGetPr(px);
                pValPoles = (double*)mxGetPr(pPoles);
                
                largox = mxGetN(pPoles)*mxGetM(pPoles);
                largo_y = largox-prepad-pospad;
                
                y_aux = mxCreateNumericMatrix(mxGetM(pPoles), mxGetN(pPoles),mxDOUBLE_CLASS, mxREAL);
                plhs[0] = mxCreateNumericMatrix(1,largo_y,mxDOUBLE_CLASS, mxREAL);
                
                yr = (double*)mxGetPr(y_aux);
                yr[prepad] = pValXr[prepad];
                
             /*Pre-pad filtering*/
                for(i = prepad-1; i>=0 ; i--){
                  
                    yr[i] = (pValXr[i+1] + pValXr[i]) + yr[i+1]*( aux_pole = pValPoles[prepad-1-i]) ;
                }
             /* Move the pointer to have indexes aligned to output indexes */
                pValPoles = pValPoles + prepad;
                
             /* forward filtering */
                for(i = 1; i<largox-prepad ; i++){
                    yr[i] = (pValXr[i-1] + pValXr[i]) + yr[i-1]* ( aux_pole = pValPoles[i]);
                    yr[i-1] = yr[i-1] + yr[i];
                }
             /* backward filtering */
                for(i = largox-prepad-2; i>=0 ; i--){
                    yr[i] += yr[i+1]*( aux_pole = pValPoles[i]);
                }
                yout_r = (double*)mxGetPr(plhs[0]);
                
            /* Gain correction  & output storage */
                for(i = 0; i<largo_y ; i++){
                    aux_pole = (1.0-(*(pValPoles++)))/2.0;
                    *(yout_r++) = (*(yr++))*aux_pole*aux_pole;/* /((double)1024)/2.0; */
                }
                
            }
            mxDestroyArray(y_aux);
            
        }
        
    }else{
        printf("iir_cqt_c\nIncorrect input/output parameters.\n\n");
        printf("See help for details\n");        
    }
    return;
}
