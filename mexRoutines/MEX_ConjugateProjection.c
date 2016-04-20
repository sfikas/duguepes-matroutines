//  res = conjugateProjection(x, weights)
// 
//  Project Kx1 vector weights .* x to hyperplane 
//  defined by weights'*y = 1, subject 
//  to the constraint sum(y) = 1, y >= 0.
//  The objective function is supposed to be of the form
//        x'Ax + bx + const.
//  A is diagonal, and 'weights' are its elements.
//  
//  Notes: 1. If weights = ones(1,K), then the result should
//            coincide to K.Blekas' BIDProjection(x).
//  Notes: 2. In the CVPR08 submission, the objective function was
//            in fact of the form x'Ax + bx + clogx + const.
//            In BLEKAS05, the objective function function was
//            x'x + bx + clogx + const.
//  G.Sfikas 21 April 2008
#define SQR(x) ((x)*(x))

#include "mex.h"
#include <stdlib.h>
#include <math.h>

#define MAXK 100
#define MALAKA printf("HOHO\n");
double K;

void ConjProjection(double *result, double *a, double *weights);
void findMaximum(double *maxValue, int *index, double x[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *x;
    double *weights;
    double *result;
    double *sizeX;

	/*/ Input*/
	x = mxGetPr(prhs[0]);
    weights = mxGetPr(prhs[1]);
    sizeX = mxGetPr(prhs[2]);
    if(x == NULL || sizeX == NULL) {
        printf("ConjProjection: Syntax is MEX_ConjuageProjection(vector, weights, number of elements in vector)\n");
        printf("ConjProjection: Example: MEX_ConjugateProjection([1 1.5 -2.1], [1 1 1], 3)\n");
    }
    K = *sizeX;
   	/*/ Output*/
	plhs[0] = mxCreateDoubleMatrix(1, K, mxREAL);
	result = mxGetPr(plhs[0]);
	ConjProjection(result, x, weights);
}


void ConjProjection(double *result, double *a, double *weights)
{
    double g[MAXK];
    double yInit[MAXK];
    double y[MAXK];
    int activeSet[MAXK];
    double sumW;
    double sumA;
    double sumL;
    int j;
    int negCount;
    double sumGINACT, sumGACTY;
    int iterations;
    double lambda[MAXK];
    double ggM; int ggExcl;
    double sumAbsX;
    
    /* Check if vector has too large variates */
    findMaximum(&ggM, &ggExcl, a);
    for(sumAbsX = 0,j = 0; j < K; j++)
        sumAbsX += fabs(a[j]);
    if(ggM > 0 && sumAbsX - ggM < 1e-5) {
        for(j = 0; j < K; j++)
            result[j] = 0;
        result[ggExcl] = 1;
        return;
    }

    for(sumW = 0, sumA = 0, j = 0; j < K; j++) {
        sumW += 1/weights[j];
        sumA += a[j];
    }
    for(j = 0; j < K; j++) 
        g[j] = -1/(sumW*weights[j]);
    //Initial projection to weights'y = 1.
    for(j = 0; j < K; j++) {
        y[j] = a[j] - g[j] + g[j]*sumA;
        yInit[j] = y[j];
    }
    // activeSet = zeros(K, 1);
    for(j = 0; j < K; j++)
        activeSet[j] = 0;
    //////////////
    // main loop
    //////////////
    for(iterations = 0; iterations < K; j++) {
        for(negCount = 0,j = 0; j < K; j++)
            if(y[j] < 0) {
                negCount++;
                activeSet[j] = 1;
            }
        if(negCount == 0) {
            for(j = 0; j < K; j++) 
                result[j] = y[j];
            return;
        }
        for(sumGACTY = 0, sumGINACT = 0, j = 0; j < K; j++) {
            if(activeSet[j] == 1)
                sumGACTY += g[j]*yInit[j];
            else
                sumGINACT += g[j];
            lambda[j] = 0;
        }
        // Iterate between recalculating Lagrange operators & proj to sum(y) = 1
        for(sumL = 0, j = 0; j < K; j++) {
            if(activeSet[j] == 0)
                lambda[j] = 0;
            else
                lambda[j] = -yInit[j] - sumGACTY / sumGINACT;
            sumL += lambda[j];
        }
        for(j = 0; j < K; j++) {
            if(activeSet[j] == 0)
                y[j] = a[j] - g[j] + g[j]*sumA + g[j]*sumL + lambda[j];
            else
                y[j] = 0;
        }
    }
    printf("conjugateProjection: Algorithm failed after K = %f iterations!\n", K);
}


void findMaximum(double *maxValue, int *index, double x[]) {
    int i;
    
    *maxValue = x[0]; *index = 0;
    for(i = 0; i < K; i++)
        if(x[i] > *maxValue) {
            *maxValue = x[i];
            *index = i;
        }
    return;
}
