//function Y = mex_matrix_op_repmat_vector(A,X,op,mode);
// Timothee Cour, 29-Aug-2006 07:49:15

/*computes Y=A op repmat(X,size(A,1),1);
 *or Y=A op repmat(X,1,size(A,2));
 *example: Y = mex_matrix_op_repmat_vector(A,X,'.*','col');
 *A:mxn
 *X:1xn (or nx1)
 *Y:mxn
 *mode:'col','row'
 */

#include <math.h>
#include <mex.h>
#include "mex_util.cpp"
#include "mex_math.cpp"


//double
void mult_double(double *A,double *X,double *Y,bool iscol,int K,int n) {
    int i,j,k=0;
    if(iscol){
        for(j=0;j<K;j++)
            for(i=0;i<n;i++)
                Y[k]=X[i]*A[k++];
    }else{
        double x;
        for(j=0;j<n;j++) {
            x=X[j];
            for(i=0;i<K;i++)
                Y[k]=x*A[k++];
        }
    }
}
void mult_double_cmplx(double *A,double *Ai,double *X,double *Xi,double *Y,double *Yi,bool iscol,int K,int n) {
    int i,j,k=0;
    if(iscol){
        for(j=0;j<K;j++)
            for(i=0;i<n;i++,k++){
                Y[k]=X[i]*A[k]-Xi[i]*Ai[k];
                Yi[k]=X[i]*Ai[k]+Xi[i]*A[k];
            }
    }else{
        double x,xi;
        for(j=0;j<n;j++) {
            x=X[j];
            xi=Xi[j];
            for(i=0;i<K;i++,k++){
                Y[k]=x*A[k]-xi*Ai[k];
                Yi[k]=x*Ai[k]+xi*A[k];
            }
        }
    }
}

void mult_self_double(double *A,double *X,bool iscol,int K,int n) {
    int i,j,k=0;
    if(iscol){
        for(j=0;j<K;j++)
            for(i=0;i<n;i++)
                A[k++]*=X[i];
    }else{
        double x;
        for(j=0;j<n;j++) {
            x=X[j];
            for(i=0;i<K;i++)
                A[k++]*=x;
        }
    }
}
void mult_self_double_cmplx(double *A,double *Ai,double *X,double *Xi,bool iscol,int K,int n) {
    int i,j,k=0;
    double temp=0;
    if(iscol){
        for(j=0;j<K;j++)
            for(i=0;i<n;i++,k++){
                //A[k++]*=X[i];
                temp=A[k];
                A[k]=A[k]*X[i]-Ai[k]*Xi[i];
                Ai[k]=temp*Xi[i]+Ai[k]*X[i];
            }
    }else{
        double x,xi;
        double temp;
        for(j=0;j<n;j++) {
            x=X[j];
            xi=Xi[j];
            for(i=0;i<K;i++,k++){
                temp=A[k];
                A[k]=A[k]*x-Ai[k]*xi;
                Ai[k]=temp*xi+Ai[k]*x;
            }
        }
    }
}

//float
void mult_float(float *A,float *X,float *Y,bool iscol,int K,int n) {
    int i,j,k=0;
    if(iscol){
        for(j=0;j<K;j++)
            for(i=0;i<n;i++)
                Y[k]=X[i]*A[k++];
    }else{
        float x;
        for(j=0;j<n;j++) {
            x=X[j];
            for(i=0;i<K;i++)
                Y[k]=x*A[k++];
        }
    }
}
void mult_float_cmplx(float *A,float *Ai,float *X,float *Xi,float *Y,float *Yi,bool iscol,int K,int n) {
    int i,j,k=0;
    if(iscol){
        for(j=0;j<K;j++)
            for(i=0;i<n;i++,k++){
                Y[k]=X[i]*A[k]-Xi[i]*Ai[k];
                Yi[k]=X[i]*Ai[k]+Xi[i]*A[k];
            }
    }else{
        float x,xi;
        for(j=0;j<n;j++) {
            x=X[j];
            xi=Xi[j];
            for(i=0;i<K;i++,k++){
                Y[k]=x*A[k]-xi*Ai[k];
                Yi[k]=x*Ai[k]+xi*A[k];
            }
        }
    }
}

void mult_self_float(float *A,float *X,bool iscol,int K,int n) {
    int i,j,k=0;
    if(iscol){
        for(j=0;j<K;j++)
            for(i=0;i<n;i++)
                A[k++]*=X[i];
    }else{
        float x;
        for(j=0;j<n;j++) {
            x=X[j];
            for(i=0;i<K;i++)
                A[k++]*=x;
        }
    }
}
void mult_self_float_cmplx(float *A,float *Ai,float *X,float *Xi,bool iscol,int K,int n) {
    int i,j,k=0;
    float temp=0;
    if(iscol){
        for(j=0;j<K;j++)
            for(i=0;i<n;i++,k++){
                //A[k++]*=X[i];
                temp=A[k];
                A[k]=A[k]*X[i]-Ai[k]*Xi[i];
                Ai[k]=temp*Xi[i]+Ai[k]*X[i];
            }
    }else{
        float x,xi;
        float temp;
        for(j=0;j<n;j++) {
            x=X[j];
            xi=Xi[j];
            for(i=0;i<K;i++,k++){
                temp=A[k];
                A[k]=A[k]*x-Ai[k]*xi;
                Ai[k]=temp*xi+Ai[k]*x;
            }
        }
    }
}

void mexFunction(int nargout, mxArray *out[], int nargin, const mxArray *in[]) {
    const mxArray *A=in[0];
    const mxArray *X=in[1];
    int n = mxGetNumberOfElements(X);
    int K = mxGetNumberOfElements(A)/n;
    mxClassID classID=mxGetClassID(A);
    bool isCmplx=mxIsComplex(A);
    mxComplexity cmplxTag;
    if (isCmplx)
        cmplxTag= mxCOMPLEX;
    else
        cmplxTag=mxREAL;
    
    char *operationStr=mxArrayToString(in[2]);
    char *modeStr=mxArrayToString(in[3]);
    
    bool iscol;
    if(strcmp(modeStr,"row")==0)
        iscol=false;
    else if(strcmp(modeStr,"col")==0)
        iscol=true;
    else
        mexErrMsgTxt("enter col or row\n");
    
    int operation=0;
    if(strcmp(operationStr,".*")==0)
        operation=0;
    else if(strcmp(operationStr,"+")==0)
        operation=1;
    else
        mexErrMsgTxt("enter +,.*\n");
    
    if(operation!=0)
        mexErrMsgTxt("todo\n");
    
    
    bool isInPlace=true;
    if(nargout==1)
        isInPlace=false;
    
    if(!isInPlace) {
        //out[0]=mxCreateDoubleMatrix(K,n,mxREAL);
        out[0]=mxCreateNumericArray(mxGetNumberOfDimensions(A), mxGetDimensions(A), classID, cmplxTag);
        mxArray *Y=out[0];
        switch(classID) {
            case mxDOUBLE_CLASS:
                if(isCmplx)
                    mult_double_cmplx(mxGetPr(A),mxGetPi(A),mxGetPr(X),mxGetPi(X),mxGetPr(Y),mxGetPi(Y),iscol,K,n);
                else
                    mult_double(mxGetPr(A),mxGetPr(X),mxGetPr(Y),iscol,K,n);
                break;
            case mxSINGLE_CLASS:
                if(isCmplx)
                    mult_float_cmplx((float*)mxGetData(A),(float*)mxGetImagData(A),(float*)mxGetData(X),(float*)mxGetImagData(X),(float*)mxGetData(Y),(float*)mxGetImagData(Y),iscol,K,n);
                else
                    mult_float((float*)mxGetData(A),(float*)mxGetData(X),(float*)mxGetData(Y),iscol,K,n);
                break;
                default:
                    mexErrMsgTxt("wrong input class type\n");
        }
    } else {
        switch(classID) {
            case mxDOUBLE_CLASS:
                if(isCmplx)
                    mult_self_double_cmplx(mxGetPr(A),mxGetPi(A),mxGetPr(X),mxGetPi(X),iscol,K,n);
                else
                    mult_self_double(mxGetPr(A),mxGetPr(X),iscol,K,n);
                break;
            case mxSINGLE_CLASS:
                if(isCmplx)
                    mult_self_float_cmplx((float*)mxGetData(A),(float*)mxGetImagData(A),(float*)mxGetData(X),(float*)mxGetImagData(X),iscol,K,n);
                else
                    mult_self_float((float*)mxGetData(A),(float*)mxGetData(X),iscol,K,n);
                break;
                default:
                    mexErrMsgTxt("wrong input class type\n");
        }
    }
}
