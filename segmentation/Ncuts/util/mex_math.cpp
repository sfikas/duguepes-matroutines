/*================================================================
// Timothee Cour, 29-Aug-2006 07:49:15

mex_math = used by a couple of mex functions 
 *=================================================================*/
 
# include "math.h"

int round2(double x) {
    //return floor(x+0.5);
    return x>=0 ? (int)(x+0.5) : (int) (x-0.5);
}
/* Problem: when compiling on opteron, says error: new declaration
int round(double x) {
    //return floor(x+0.5);
    return x>=0 ? (int)(x+0.5) : (int) (x-0.5);
}*/

double min(double x,double y) {
    return x<y?x:y;
}
double max(double x,double y) {
    return x>y?x:y;
}
int max(int x,int y) {
    return x>y?x:y;
}
int min(int x,int y) {
    return x<y?x:y;
}

double vec_max(double *x,int n) {
    double res=x[0];
    for(int i=1;i<n;i++)
        if(res<x[i])
            res=x[i];
    return res;
}
int vec_max(int *x,int n) {
    int res=x[0];
    for(int i=1;i<n;i++)
        if(res<x[i])
            res=x[i];
    return res;
}
int vec_min(int *x,int n) {
    int res=x[0];
    for(int i=1;i<n;i++)
        if(res>x[i])
            res=x[i];
    return res;
}
double vec_min(double *x,int n) {
    double res=x[0];
    for(int i=1;i<n;i++)
        if(res>x[i])
            res=x[i];
    return res;
}
double dot(double *x,double *y, int n) {
    double temp = 0;
    for(int i=0;i!=n;i++)
        temp+=*x++ * *y++;
    return temp;
}

void scalar_times_vec(double a, double *x,double *y, int n) {
    // y=a*x
    for(int i=0;i!=n;i++)
        *y++ = *x++ * a;
}
void scalar_plus_vec_self(int a, int *x, int n) {
    // x=a+x
    for(int i=0;i!=n;i++)
        *x++ += a;
}

void vec_plus_vec(double *x1, double *x2,double *y, int n) {
    // y=x1+x2
    for(int i=0;i!=n;i++)
        y[i] = x1[i] + x2[i];
        //*y++ = *x1++ + *x2++;
}


void ind2sub(int*ind,int*indi,int*indj,int p,int q,int n){
	//caution: output is in matlab conventions//TODO;correct this [mex_ind2sub
    int i;
    int*pindi=indi;
    int*pindj=indj;
    int*pind=ind;
    int temp;
    for(i=0;i<n;i++){
        temp=*pind++-1;
        *pindi++=temp%p+1;
        *pindj++=temp/p+1;
    }
}

//void symmetrize_sparse()


