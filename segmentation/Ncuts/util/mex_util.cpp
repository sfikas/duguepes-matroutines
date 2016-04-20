/*================================================================
// Timothee Cour, 29-Aug-2006 07:49:15

mex_util = used by a couple of mex functions 
*=================================================================*/

//#undef NDEBUG

#include <sstream>
# include "math.h"
//#include <string>
#include <string.h>
#include <time.h>
#include "mat.h"
#include <vector>
#include <map>
#include <string>
#include <cstdarg>
//#include <assert.h>

using namespace std;

//MACROS


//# define assert(isOK) ( (isOK) ? (void)0 : (void)mexErrMsgIdAndTxt("","!! Assert '%s' failed on line %d in file '%s'\n", #isOK, __LINE__, __FILE__) )

#define printLine mexPrintf("line %d in file '%s'\n", __LINE__, __FILE__)
#define printError mexErrMsgIdAndTxt("","line %d in file '%s'\n", __LINE__, __FILE__)

#ifndef NDEBUG
# define ASSERT( isOK ) \
	( (isOK) ? \
	(void)0 : \
	(void)mexErrMsgIdAndTxt("","!! Assert '%s' failed on line %d in file '%s'\n", #isOK, __LINE__, __FILE__) )
#else
# define ASSERT( unused ) do {} while( false )
#endif

# define assert( isOK ) \
	( (isOK) ? \
	(void)0 : \
	(void)mexErrMsgIdAndTxt("","!! Assert '%s' failed on line %d in file '%s'\n", #isOK, __LINE__, __FILE__) )

void disp(string s1,string s2);

#define disp2(x)(disp(x,#x))

#ifndef NDEBUG
# define DEBUG(val){mexPrintf("line %d, file '%s': ", __LINE__, __FILE__);disp2(val);}
#else
//# define DEBUG(unused) {do {} while( false )}
# define DEBUG(unused) {}
#endif


void getSizes(const mxArray *A,int &p,int &q,int &k){
	if (mxGetNumberOfDimensions(A)<=2){
		p=mxGetM(A);
		q=mxGetN(A);
		k=1;
	}
	else {
		const int *dims = mxGetDimensions(A);
		p=dims[0];
		q=dims[1];
		k=dims[2];
	}
}
int getSize3(const mxArray *A){
	int k=0;
	if (mxGetNumberOfDimensions(A)<=2)
		k = 1; 
	else {
		const int *dims = mxGetDimensions(A);
		k = dims[2];
	}
	return k;
}

void** deal_output(int num,...) {
	//deprecated:use deal_outputs
	void**output=(void**)mxCalloc(num,sizeof(void*));
	va_list arguments;                     
	va_start(arguments,num);           
	for(int i=0;i<num;i++)        
		output[i]=va_arg(arguments,void*); 
	va_end(arguments);                  
	return output;
}
vector<void*> deal_outputs(int num,...) {
	vector<void*>output(num);
	//void**output=(void**)mxCalloc(num,sizeof(void*));
	va_list arguments;                     
	va_start(arguments,num);           
	for(int i=0;i<num;i++)        
		output[i]=va_arg(arguments,void*); 
	va_end(arguments);                  
	return output;
}
mxArray* createMxArrayInt(int m,int n){
	int ndim=2;
	const int dims[]={m,n};
	return mxCreateNumericArray(ndim,dims,mxINT32_CLASS,mxREAL);
}

void fill_out_full_int(mxArray *out[],int nout,int *val,int m,int n) {
	out[nout] = mxCreateDoubleMatrix(m, n, mxREAL);
	double *temp = mxGetPr(out[nout]);
	int nnz = m*n;
	for (int i=0;i<nnz;i++)
		temp[i] = (double)val[i];
}
void fill_out_full(mxArray *out[],int nout,double *val,int m,int n) {
	out[nout] = mxCreateDoubleMatrix(m, n, mxREAL);
	double *temp = mxGetPr(out[nout]);
	int nnz = m*n;
	for (int i=0;i<nnz;i++)
		temp[i] = val[i];
}

void readSparse(const mxArray *A,double**pr,int**ir,int**jc,int *m,int *n){
	*pr=mxGetPr(A);
	*ir=mxGetIr(A);
	*jc=mxGetJc(A);
	*m=mxGetM(A);
	*n=mxGetN(A);
}

int* readInt32(const mxArray *A) {
	switch (mxGetClassID(A))  {
		case mxINT32_CLASS: return (int*)mxGetData(A);  break;//TODO:attention with subsequent modif of pointer
		case mxDOUBLE_CLASS: {
			double *pr = mxGetPr(A);
			int n = mxGetM(A)*mxGetN(A);
			int *res;
			if(n>0)
				res = (int*)mxCalloc(n,sizeof(int));
			else
				res = NULL;
			for (int i=0;i<n;i++)
				res[i] = (int)pr[i];
			return res;
							 }
		default: printError;//mexErrMsgTxt("bad input");
			return NULL;
	}
}

vector<string>getFieldNames(const mxArray*A){
	int n=mxGetNumberOfFields(A);
	vector<string> fields;
	fields.resize(n);
	for(int i=0;i<n;i++)
		fields[i]=mxGetFieldNameByNumber(A,i);
	return fields;
}
bool isField(const mxArray*A,string field){
	vector<string>fields=getFieldNames(A);
	for(int i=0;i<fields.size();i++)
		if(fields[i]==field)
			return true;
	return false;
}
//mxArray*getFieldByName(const mxArray*A,int i, const char *fname){
mxArray*getFieldByName(const mxArray*A,int i, string fname){
	assert(isField(A,fname));//inefficient?
	int field=mxGetFieldNumber(A,fname.c_str());
	return mxGetFieldByNumber(A,i,field);
}
//void setFieldByName(mxArray*A,int i, const char *fname,mxArray*B){
void setFieldByName(mxArray*A,int i, string fname,mxArray*B){
	assert(isField(A,fname));//inefficient?
	int field=mxGetFieldNumber(A,fname.c_str());
	mxSetFieldByNumber(A,i,field,B);      
}
mxArray* array2mxArray(double*x,int m,int n) {
	mxArray* A=mxCreateDoubleMatrix(m, n, mxREAL);
	double *pA=mxGetPr(A);
	int mn=m*n;    
	for(int i=0;i<mn;i++,pA++,x++)
		*pA=*x;
	x-=mn;
	return A;
}
mxArray* array2mxArray(bool*x,int m,int n) {
	mxArray* A=mxCreateLogicalMatrix(m, n);
	bool *pA=(bool*)mxGetData(A);
	int mn=m*n;    
	for(int i=0;i<mn;i++,pA++,x++)
		*pA=*x;
	x-=mn;
	return A;
}
mxArray* array2mxArray(int*x,int m,int n) {
	int ndim=2;
	const int dims[]={m,n};
	mxArray* A=mxCreateNumericArray(ndim,dims,mxINT32_CLASS,mxREAL);
	int *pA=(int*)mxGetData(A);
	int mn=m*n;    
	for(int i=0;i<mn;i++,pA++,x++)
		*pA=*x;
	x-=mn;
	return A;
}
mxArray* create_mxArray(int m,int n,mxClassID classID) {
	int ndim=2;
	const int dims[]={m,n};
	mxArray* A=mxCreateNumericArray(ndim,dims,classID,mxREAL);
	return A;
}
mxArray* array2mxArray(int*x,int m,int n,mxClassID classID) {
	int ndim=2;
	const int dims[]={m,n};
	mxArray* A=mxCreateNumericArray(ndim,dims,classID,mxREAL);
	int mn=m*n;    
	switch(classID){
		case mxDOUBLE_CLASS:{
			double *pA=mxGetPr(A);
			for(int i=0;i<mn;i++,pA++,x++)
				*pA=(double)*x;
			x-=mn;
			break;
							}
		default:
			mexErrMsgTxt("wrong argument\n");
	}
	return A;
}
mxArray* array2mxArray(double*x,int m,int n,int k,mxClassID classID) {
	int ndim=3;
	const int dims[]={m,n,k};
	mxArray* A=mxCreateNumericArray(ndim,dims,classID,mxREAL);
	int mnk=m*n*k;    
	switch(classID){
		case mxDOUBLE_CLASS:{
			double *pA=mxGetPr(A);
			for(int i=0;i<mnk;i++,pA++,x++)
				*pA=(double)*x;
			x-=mnk;
			break;
							}
		default:
			mexErrMsgTxt("wrong argument\n");
	}
	return A;
}


mxArray* array2mxArray(vector<double> &v){
	int n=v.size();
	mxArray* A;
	//if(n==0)
	//	A=mxCreateDoubleMatrix(0, 0, mxREAL);
	//else
	//	A=mxCreateDoubleMatrix(n, 1, mxREAL);

	A=mxCreateDoubleMatrix(n, 1, mxREAL);
	double *pA=mxGetPr(A);
	for(int i=0;i<n;i++)
		*pA++=v[i];
	return A;
}
mxArray* array2mxArray(vector<bool> &v){
	int n=v.size();
	mxArray* A=mxCreateLogicalMatrix(n, 1);
	bool *pA=(bool*)mxGetData(A);
	for(int i=0;i<n;i++)
		*pA++=v[i];
	return A;
}
mxArray* array2mxArray(vector<int> &v,mxClassID classID,int a){
	int n=v.size();
	vector<int>::iterator pv=v.begin();
	switch (classID)  {
		case mxINT32_CLASS:			
			{
				mxArray* A=createMxArrayInt(n,1);
				int*pA=(int*)mxGetData(A);
				for(int i=0;i<n;i++)
					*pA++=a+*pv++;
				return A;  
				break;
			}						
		case mxDOUBLE_CLASS: 
			{
				double a2=(double)a;
				mxArray* A=mxCreateDoubleMatrix(n,1,mxREAL);
				double*pA=mxGetPr(A);
				for(int i=0;i<n;i++)
					*pA++=a2+(double)*pv++;
				return A;  
				break;										
			}			
		default: printError;//mexErrMsgTxt("bad input");
			return NULL;
	}
}
void mxArray2array(const mxArray* A,vector<int> &v,int a){
	int n=mxGetNumberOfElements(A);
	v.resize(n);
	vector<int>::iterator pv=v.begin();
	switch (mxGetClassID(A)) {
		case mxINT32_CLASS: 
			{
				int*pA=(int*)mxGetData(A);
				for(int i=0;i<n;i++)
					*pv++=(a+*pA++);			
				break;							
			}
		case mxDOUBLE_CLASS: 
			{
				double*pA=mxGetPr(A);
				for(int i=0;i<n;i++)
					*pv++=(int)(a+*pA++);
				break;							
			}
		default: printError;//mexErrMsgTxt("bad input");
	}	
}
void mxArray2array(const mxArray* A,vector<double> &v,double a){
	int n=mxGetNumberOfElements(A);
	v.resize(n);
	vector<double>::iterator pv=v.begin();
	switch (mxGetClassID(A)) {
		case mxINT32_CLASS: 
			{
				int*pA=(int*)mxGetData(A);
				for(int i=0;i<n;i++)
					*pv++=(double)(a+*pA++);			
				break;							
			}
		case mxDOUBLE_CLASS: 
			{
				double*pA=mxGetPr(A);
				for(int i=0;i<n;i++)
					*pv++=(double)(a+*pA++);
				break;							
			}
		default: printError;//mexErrMsgTxt("bad input");
	}	
}
void mxArray2array(const mxArray* A,vector<bool> &v){
	int n=mxGetNumberOfElements(A);
	v.resize(n);
	vector<bool>::iterator pv=v.begin();
	switch (mxGetClassID(A)) {
		case mxLOGICAL_CLASS:
			{
				bool*pA=(bool*)mxGetData(A);
				for(int i=0;i<n;i++)
					*pv++=*pA++;	
				break;							
			}
		default: printError;
	}	
}
void mxArray2array(const mxArray* A,string &v){
	assert(mxGetClassID(A)==mxCHAR_CLASS);
	char *v0=mxArrayToString(A);
	v=v0;
}
void mxArrayStruct2array(const mxArray* A,int i,string field,string &v){
	const mxArray*Ai=getFieldByName(A,i,field);
	mxArray2array(Ai,v);
}
void mxArrayStruct2array(const mxArray* A,int i,string field,vector<double> &v,double a){
	const mxArray*Ai=getFieldByName(A,i,field);
	mxArray2array(Ai,v,a);
}
void mxArrayStruct2array(const mxArray* A,int i,string field,vector<int> &v,int a){
	const mxArray*Ai=getFieldByName(A,i,field);
	mxArray2array(Ai,v,a);
}
void mxArrayStruct2array(const mxArray* A,int i,string field,int &v,int a){
	const mxArray*Ai=getFieldByName(A,i,field);
	v=a+(int)mxGetScalar(Ai);
}
void mxArrayStruct2array(const mxArray* A,int i,string field,double &v,double a){
	const mxArray*Ai=getFieldByName(A,i,field);
	v=a+(double)mxGetScalar(Ai);
}


double* createArray(int n,double val){
	if(n==0)
		return NULL;
	double *x = (double*)mxCalloc(n,sizeof(double));
	for(int i=0;i!=n;i++,x++)
		*x=val;
	x-=n;
	return x;
}

bool* createArray(int n,bool val){
	if(n==0)
		return NULL;
	bool *x = (bool*)mxCalloc(n,sizeof(bool));
	for(int i=0;i!=n;i++,x++)
		*x=val;
	x-=n;
	return x;
}
int* createArray(int n,int val){
	if(n==0)
		return NULL;
	int *x = (int*)mxCalloc(n,sizeof(int));
	for(int i=0;i!=n;i++,x++)
		*x=val;
	x-=n;
	return x;
}

mxArray*createStructArray(int n, vector<string>fields){
	int nFields=fields.size();
	const char **field_names = (const char **)mxCalloc(nFields, sizeof(*field_names));
	for(int i=0;i<nFields;i++)
		field_names[i]=fields[i].c_str();
	mxArray*A=mxCreateStructMatrix(n,1,nFields,field_names);
	mxFree(field_names);
	return A;
}
//bool isEmpty(const mxArray*A){
//	return mxGetNumberOfElements(A)==0;
//}




/* use mxArrayToString
char *str readString(const mxArray *A) {
int buflen = (mxGetM(A) * mxGetN(A)) + 1;
char *str=mxCalloc(buflen, sizeof(char));
int status = mxGetString(A, str, buflen);
return str;
}
*/

void copyDims(const mxArray *A,mxArray *B){
	int number_dims= mxGetNumberOfDimensions(A);
	const int* dims=mxGetDimensions(A);
	mxSetDimensions(B,dims, number_dims);
}
void setDimensions(mxArray *A,int m,int n){
	int ndim=2;
	const int dims[]={m,n};
	int isError = mxSetDimensions(A, dims, ndim);
	if(isError)
		printError;
}

void setDimensions(mxArray *A,int m,int n,int k){
	int ndim=3;
	const int dims[]={m,n,k};
	int isError = mxSetDimensions(A, dims, ndim);
	if(isError)
		printError;
}

void swap(double *x,double *y) {
	double temp=*x;
	*x=*y;
	*y=temp;
}

double *copyArray(const mxArray *A) {
	//TODO:use mxDuplicateArray?
	int n=mxGetM(A)*mxGetN(A);
	return (double*) memcpy(mxCalloc(n, sizeof(double)), mxGetPr(A), sizeof(double)*n);
}
double *copyArray(double *x,int n) {
	return (double*) memcpy(mxCalloc(n, sizeof(double)), x, sizeof(double)*n);
}
void copyArray(double *x,int n,double *y) {
	memcpy(y, x, sizeof(double)*n);
}


//assumes indi,indj ordered
void compute_j2Jc(int*jc,int*indj,int n,int nnz) {    
	int k = 0;    
	for(int j=0;j<n;j++) {
		jc[j] = k;
		while(k<nnz && indj[k]==j)
			k++;
	}
	jc[n] = k;
}
/*
int compute_ij2k(int*jc,int*indj,int n,int nnz) {
int k = 0;    
for(int j=0;j<n;j++) {
jc[j] = k;
while(k<nnz && indj[k]==j) {
k++;
}        
}
jc[n] = k;
}
*/

mxArray* sparseSorted(int*indi,int*indj,double*val,int m,int n,int nnz){
	mxArray*W=mxCreateSparse(m,n,nnz,mxREAL);
	double *pr = mxGetPr(W);
	int *ir = mxGetIr(W);
	int *jc = mxGetJc(W);
	int k;
	int*pir=ir,*pindi=indi;
	double*ppr=pr,*pval=val;
	for(k=0;k<nnz;k++){
		*pir++=*pindi++;
		*ppr++=*pval++;
	}
	compute_j2Jc(jc,indj,n,nnz);
	return W;
}
mxArray* sparseSorted(int*indi,int*indj,double val,int m,int n,int nnz){
	mxArray*W=mxCreateSparse(m,n,nnz,mxREAL);
	double *pr = mxGetPr(W);
	int *ir = mxGetIr(W);
	int *jc = mxGetJc(W);
	int k;
	int*pir=ir,*pindi=indi;
	double*ppr=pr;
	for(k=0;k<nnz;k++){
		*pir++=*pindi++;
		*ppr++=val;
	}
	compute_j2Jc(jc,indj,n,nnz);
	return W;
}

/*
void fill_out_sparse(mxArray *out[],int nout,int *pr, int *jc,double *pr,int nnz, int m,int n) {
out[nout] = mxCreateSparse(m, n, nnz, mxREAL);
double *temp = mxGetPr(out[nout]);
for (int i=0;i<nnz;i++)
temp[i] = val[i];
}*/


void displayMatrix(double *W,int m,int n) {
	mxArray *W2 = mxCreateDoubleMatrix(m, n, mxREAL);
	double *prW = mxGetPr(W2);
	int k;
	for (k=0;k<m*n;k++)
		prW[k] = W[k];
	mxArray *args[1];
	args[0] = (mxArray*) W2;
	mexCallMATLAB(0, NULL, 1, args, "disp");
	mxDestroyArray(W2);
}
void displayMatrix(int *W,int m,int n) {
	mxArray *W2 = mxCreateDoubleMatrix(m, n, mxREAL);
	double *prW = mxGetPr(W2);
	int k;
	for (k=0;k<m*n;k++)
		prW[k] = (double)(W[k]);
	mxArray *args[1];
	args[0] = (mxArray*) W2;
	//NOTE: syntax will not work with -argcheck compiler option
	mexCallMATLAB(0, NULL, 1, args, "disp");
	mxDestroyArray(W2);
}

void imagescMatrix(double *W,int m,int n) {
	mxArray *W2 = mxCreateDoubleMatrix(m, n, mxREAL);
	double *prW = mxGetPr(W2);
	int k;
	for (k=0;k<m*n;k++)
		prW[k] = W[k];
	mxArray *args[1];
	args[0] = (mxArray*) W2;
	mexCallMATLAB(0, NULL, 0, NULL, "figure2");
	mexCallMATLAB(0, NULL, 1, args, "imagesc");
	mxDestroyArray(W2);
}
void showImages(double *W,int m,int n,int k) {
	int dims[]={m,n,k};
	mxArray *W2 = mxCreateNumericArray(3,dims,  mxDOUBLE_CLASS,mxREAL);
	double *prW = mxGetPr(W2);
	int i;
	for (i=0;i<m*n*k;i++)
		prW[i] = W[i];
	mxArray *args[1];
	args[0] = (mxArray*) W2;
	mexCallMATLAB(0, NULL, 0, NULL, "figure2");
	mexCallMATLAB(0, NULL, 1, args, "showImages");
	mxDestroyArray(W2);
}

void save2matfile(mxArray *A,const char *file,const char *varname){
	MATFile *pmat = matOpen(file, "w");
	if (pmat == NULL)
		mexErrMsgTxt("error accessing mat file\n");
	int status = matPutVariable(pmat, varname, A);
	if (status != 0)
		mexErrMsgTxt("error writing mat file\n");
	if (matClose(pmat) != 0)
		printf("Error closing file %s\n",file);      
}

void debug(bool condition,int val) {
	if(condition) {
		mexPrintf("val1 = %d\n",val);
	}
}
void debug(bool condition,int val,char *s) {
	if(condition) {
		//mexPrintf(s+" = %d\n",val);
		//mexPrintf(strcat(s," = %d\n"),val);
		printf(strcat(s," = %d\n"),val);
	}
}
void debug(bool condition,double val) {
	if(condition) {
		mexPrintf("val1 = %1.3g\n",val);
	}
}
void stop() {
	mexErrMsgTxt("Stop\n");
}
void stop(char *s) {
	mexPrintf("%s\n",s);    
	mexErrMsgTxt("Stop\n");
}
//void stop(int i,char *s) {
//	//stop(__LINE__,__FILE__);
//	//mexPrintf("Stopped at: File: %s, line: %d\n",__FILE__,__LINE__);
//	mexPrintf("Stopped at: File: %s, line: %d\n",s,i);
//	mexErrMsgTxt("Stop\n");
//}

void checkBounds(int val,int lim) {
	if(val<0) {		
		mexPrintf("checkBounds failed: val = %d < 0\n",val);
		mexErrMsgTxt("Error\n");
	}
	if(val>=lim){
		mexPrintf("checkBounds failed: val = %d >= %d\n",val,lim);
		mexErrMsgTxt("Error\n");
	}
}

void disp(string s,string var) {
	mexPrintf("%s=%s\n",var.c_str(),s.c_str());
}
void disp(string s) {
	mexPrintf("%s\n",s.c_str());
}
void disp(int x) {
	mexPrintf("? = %d\n",x);
}
void disp(double x) {
	mexPrintf("? = %1.3g\n",x);
}
void disp(int x,char *s) {
	mexPrintf("%s = %d\n",s,x);
	//mexPrintf(strcat(s," = %d\n"),x);
}
void disp(double x,char *s) {
	mexPrintf("%s = %1.3g\n",s,x);
	//mexPrintf(strcat(s," = %1.3g\n"),x);
}


void disp(vector<double> &x){
	disp("vector:");
	for(int i=0;i<x.size();i++)
		mexPrintf("?[%d] = %1.3g\n",i,x[i]);
}
void disp(vector<int> &x){
	disp("vector:");
	for(int i=0;i<x.size();i++)
		mexPrintf("?[%d] = %d\n",i,x[i]);
}
void disp(vector<vector<int> > &x){
	disp("vector:");
	for(int j=0;j<x.size();j++){
		mexPrintf("?[%d]:",j);
		for(int i=0;i<x[j].size();i++)
			mexPrintf("?[%d][%d] = %d ; ",j,i,x[j][i]);
		mexPrintf("\n");
	}
}

void dispTime(long time){
	mexPrintf("time = %10.5g\n",((double)time)/CLOCKS_PER_SEC);   
}
void dispTime(long time[],int n){
	double sum=0,timei;
	for(int i=0;i<n;i++){
		timei=((double)time[i])/CLOCKS_PER_SEC;
		mexPrintf("time[%d] = %10.5g\n",i,timei);   
		sum+=timei;
	}
	mexPrintf("TOTAL = %10.5g\n",sum);   
}
void dispTime(map<string,long>&timing){
	double sum=0,timei;

	for(map<string,long>::iterator iter = timing.begin(); iter != timing.end(); iter++){
		timei=((double)(*iter).second)/CLOCKS_PER_SEC;
		mexPrintf("timing[%s] = %10.5g\n",(*iter).first.c_str(),timei);   
		sum+=timei;
    }
	mexPrintf("TOTAL = %10.5g\n",sum);   
}
template<class T>void deleteVector(vector<T> &v){
	int n=v.size();
	for(int i=0;i<n;i++)
		delete v[i];
}
template<class T>void array2vector(T*x,vector<T> &v,int n){
	v.resize(n);
	vector<T>::iterator pv=v.begin();
	for(int i=0;i<n;i++)
		*pv++=*x++;
}
template<class T>void vector2array(vector<T> &v,T*x){
	int n=v.size();
	vector<T>::iterator pv=v.begin();
	for(int i=0;i<n;i++)
		*x++=*pv++;
}





string toStr(int x){
	ostringstream s0;
	s0 << x;
	string s = s0.str();
	return s;
}
string toStr(double x){
	ostringstream s0;
	s0 << x;
	string s = s0.str();
	return s;
}
string toStr(long x){
	ostringstream s0;
	s0 << x;
	string s = s0.str();
	return s;
}
