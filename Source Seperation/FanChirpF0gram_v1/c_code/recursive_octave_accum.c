/* This function accumulates the elements of a matrix columnwise that are 
 * at a distance multiple of n_octave
 *             y(i,:) = \sum_k x(i+k*n_octave,:)
 */

#include "mex.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*Declarations*/
    mxArray *px;
    long n_octave;
    double *yr;
    double* y_aux;
    int filas;
    int i;
    int j;
    if(nrhs == 2 && nlhs ==1){
        /* px is the pointer to the matrix to accumulate*/
        px = (mxArray*)prhs[0];
        
        if(mxIsClass(px, "double")){
            double *pValXr;
            
            long cols;
            
            pValXr = (double*)mxGetPr(px);
             /* prhs[1] is the distance between two values thay are one octave appart*/
            n_octave = (long)(*((double*)mxGetPr((mxArray*)prhs[1])));
            
            filas = mxGetM(px);
            cols = mxGetN(px);
            
            /*Salida*/
            plhs[0] = mxCreateNumericMatrix(mxGetM(px), mxGetN(px),mxDOUBLE_CLASS, mxREAL);
            
            yr = (double*)mxGetPr(plhs[0]);

            /*copy to operate inplace*/
            memcpy(yr,pValXr,filas*cols*8);
            
            for(j = 0; j < cols ; j++){ /* for each column */
                y_aux = yr + j*filas;
                /* for all the elements that its octave is in the matrix*/
                for(i = filas-n_octave-1; i>=0 ; i--){ 
                    y_aux[i] += y_aux[i+n_octave]; /*accumulate (recursively as it is inplace)*/
                }
            }
        }
    }else{
        printf("recursive_octave_accum.m\nIncorrect input/output parameters.\n\n");
        printf("See help. \n");
    }
    return;
}
