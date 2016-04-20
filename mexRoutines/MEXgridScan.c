/*//////////////////////
//
// MEXmultiResolution
//
// G.Sfikas 14 Mar 2008
////////////////////////*/

#define EPS 10e-10
#define DIRECTIONS 2
#define NHOODSIZE 4
#define SQR(x) ((x)*(x))
#define MAXK 100
#define PAPARA printf("PAPARA!\n");


#include "mex.h"
#include <stdlib.h>
#include <math.h>

int imageSize1;
int imageSize2;
int K;
double _mrLevel;

void ConjProjection(double *result, double *a, double *weights);
void BIDProjection(double *result, double *x);
double solveQuad(double a, double b, double c);
void findMaximum(double *maxValue, int *index, double x[]);
void mainFunc(double mrLevel, double *w, double *z, double *u, double *beta, double flagMultiResAlsoZ);
int optimizeGrid(double *result, int *gridLimits1, int *gridLimits2, double *w, double *z, double *u, double *beta,
    double flagMultiResAlsoZ);
void getNeighbourIndex(int **gridBordersIndex, int **neighboursIndex, int *noNeighbours, int n, int *gL1, int *gL2);
int getNeighbourInfo(double *pos, int n);
long acc(int dir1, int dir2, int dir3);
long accB(int dir1, int dir2);
long accU(int dir1, int dir2, int dir3);



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *mrLevel;
    double *imageSize;
    double *w;
    double *z;
    double *u;
    double *beta;
    double *nKernels;
    double *newW;
    double *flagMultiResAlsoZ;

	/*/ Input*/
	mrLevel = mxGetPr(prhs[0]);
    imageSize = mxGetPr(prhs[1]);
    w = mxGetPr(prhs[2]);
    z = mxGetPr(prhs[3]);
    u = mxGetPr(prhs[4]);
    beta = mxGetPr(prhs[5]);
    nKernels = mxGetPr(prhs[6]);
    flagMultiResAlsoZ = mxGetPr(prhs[7]);
    /* */
    imageSize1 = imageSize[0];
    imageSize2 = imageSize[1];
    K = *nKernels;
   	/*/ Output*/
	plhs[0] = mxCreateDoubleMatrix(imageSize1*imageSize2, K, mxREAL);
	newW = mxGetPr(plhs[0]);

    memcpy(newW, w, imageSize1*imageSize2*K*sizeof(double));
	mainFunc(*mrLevel, newW, z, u, beta, *flagMultiResAlsoZ);
}

void mainFunc(double mrLevel, double *w, double *z, double *u, double *beta, double flagMultiResAlsoZ) {
    double normalGridSize;
    int gridCount;
    int gridSuccessCount;
    double optimizationResult[MAXK];
    int m1, m2;
    int gridLimits1[2]; int gridLimits2[2];
    int j;
    int xx1; int xx2;
    int step;
//     double *wNew;
    
//     wNew = (double *) malloc(sizeof(double)*imageSize1*imageSize2*K);
//     memcpy(wNew, w, imageSize1*imageSize2*K);
    _mrLevel = mrLevel;
    normalGridSize = (int) pow(2.0, (double) mrLevel);
    gridCount = 0; gridSuccessCount = 0;
    step = normalGridSize/2;
//  step = normalGridSize;
    if(step == 0)
        step = 1;
    for(m1 = 1; m1 <= imageSize1 ; m1 += step)        
        for(m2 = 1; m2 <= imageSize2 ; m2 += step) {            
            gridLimits1[0] = m1;
            gridLimits2[0] = m2;
            if(m1 + normalGridSize > imageSize1)
                gridLimits1[1] = imageSize1;
            else
                gridLimits1[1] = m1 + normalGridSize - 1;
            if(m2 + normalGridSize > imageSize2)
            	gridLimits2[1] = imageSize2;
            else
            	gridLimits2[1] = m2 + normalGridSize - 1;
            if(optimizeGrid(optimizationResult, gridLimits1, gridLimits2, w, z, u, beta, flagMultiResAlsoZ) == 0)
            {
                gridSuccessCount++;                
                for(j = 1; j <= K; j++)
                    for(xx1 = gridLimits1[0]; xx1 <= gridLimits1[1]; xx1++)
                        for(xx2 = gridLimits2[0]; xx2 <= gridLimits2[1]; xx2++)
                            w[acc(xx1, xx2, j)] = optimizationResult[j-1];                            
//                             wNew[acc(xx1, xx2, j)] = optimizationResult[j-1];
            }
            gridCount++;            
   }
//    memcpy(w, wNew, imageSize1*imageSize2*K);
   printf("Pass on grid size 2to%f: %d out of %d total grids optimized\n", mrLevel, gridSuccessCount, gridCount);
}

int optimizeGrid(double *result, int *gridLimits1, int *gridLimits2, double *w, double *z, double *u, double *beta,
        double flagMultiResAlsoZ) {
    /* Return 0 for success, otherwise -1 */
    int j, k, n;
    double pos[2];
    double direction;
    double a, b, c;
    int *neighboursIndex;    /* Here'll be stored indices for current neighbours */
    int *gridBordersIndex;   /* Here'll be stored the corresponding members of the grid */
    int noNeighbours;           /* Number of neighbours */
    double energy, oldEnergy, oldProjEnergy;
    double *res;
    int r1, r2;
    double projectionWeights[MAXK];
    double oldProjResult[MAXK];
    
    res = (double *) malloc(K*sizeof(double));
    oldEnergy = energy = 0;
    for(j = 1; j<= K; j++) {
        a = 0; b = 0; 
        if(flagMultiResAlsoZ == 0)
            c = 0;
        else
            c = 1;
        for(k = 1; k <= NHOODSIZE; k++) {
            direction = getNeighbourInfo(pos, k);
            getNeighbourIndex(&gridBordersIndex, &neighboursIndex, &noNeighbours, k, gridLimits1, gridLimits2);
            for(n = 0; n < noNeighbours; n++) {
                a -= u[accU(k, j, gridBordersIndex[n])] / beta[accB(direction, j)]; 
                b += u[accU(k, j, gridBordersIndex[n])]*w[neighboursIndex[n] + imageSize1*imageSize2*(j-1)] / beta[accB(direction, j)];
                /* At the same time, compute the energy prior to the (grid) optimization */
                oldEnergy += u[accU(k, j, gridBordersIndex[n])]* 
                    SQR(w[neighboursIndex[n] + imageSize1*imageSize2*(j-1)] - w[gridBordersIndex[n] + imageSize1*imageSize2*(j-1)])
                        / beta[accB(direction, j)];
            }
            if(neighboursIndex != NULL) {
                free(neighboursIndex);
                free(gridBordersIndex);
            }
        }
        /**/
        for(r1 = gridLimits1[0]; r1 <= gridLimits1[1]; r1++)
            for(r2 = gridLimits2[0]; r2 <= gridLimits2[1]; r2++) {
                oldEnergy -= z[acc(r1, r2, j)]*log(w[acc(r1, r2, j)] + EPS);
                if(flagMultiResAlsoZ == 0)
                    c += z[acc(r1, r2, j)];
                else
                    c *= z[acc(r1, r2, j)];
            }
        c *= 0.5;
        if(flagMultiResAlsoZ == 1) //Number of grid members
            c *= (gridLimits1[1] - gridLimits1[0] + 1)*(gridLimits2[1] - gridLimits2[0] + 1);
        res[j-1] = solveQuad(a, b, c);
        // experimental:
        projectionWeights[j - 1] = -a;
    }
    ConjProjection(result, res, projectionWeights);
    /* Check if the "result" gives a lower energy/cost (higher likelihood) 
     * Return 0 if it does, otherwise -1. */
    for(j = 1; j<= K; j++) {
        for(k = 1; k <= NHOODSIZE; k++) {
            direction = getNeighbourInfo(pos, k);
            getNeighbourIndex(&gridBordersIndex, &neighboursIndex, &noNeighbours, k, gridLimits1, gridLimits2);
            for(n = 0; n < noNeighbours; n++) {
                energy += u[accU(k, j, gridBordersIndex[n])]* 
                    SQR(w[neighboursIndex[n] + imageSize1*imageSize2*(j-1)] - result[j-1]) / beta[accB(direction, j)];
            }
            if(neighboursIndex != NULL) {
                free(neighboursIndex);
                free(gridBordersIndex);
            }
        }
        /**/
        for(r1 = gridLimits1[0]; r1 <= gridLimits1[1]; r1++)
            for(r2 = gridLimits2[0]; r2 <= gridLimits2[1]; r2++)
                energy -= z[acc(r1, r2, j)]*log(result[j-1] + EPS);
    }
    /// The following section is for test purpose, and can be omitted.
    BIDProjection(oldProjResult, res);
    oldProjEnergy = 0;
    for(j = 1; j<= K; j++) {
        for(k = 1; k <= NHOODSIZE; k++) {
            direction = getNeighbourInfo(pos, k);
            getNeighbourIndex(&gridBordersIndex, &neighboursIndex, &noNeighbours, k, gridLimits1, gridLimits2);
            for(n = 0; n < noNeighbours; n++) {
                oldProjEnergy += u[accU(k, j, gridBordersIndex[n])]* 
                    SQR(w[neighboursIndex[n] + imageSize1*imageSize2*(j-1)] - oldProjResult[j-1]) / beta[accB(direction, j)];
            }
            if(neighboursIndex != NULL) {
                free(neighboursIndex);
                free(gridBordersIndex);
            }
        }
        /**/
        for(r1 = gridLimits1[0]; r1 <= gridLimits1[1]; r1++)
            for(r2 = gridLimits2[0]; r2 <= gridLimits2[1]; r2++)
                oldProjEnergy -= z[acc(r1, r2, j)]*log(oldProjResult[j-1] + EPS);
    }
    /// Until here
    if(oldProjEnergy < energy) {
        for(j = 0; j < K; j++)
            result[j] = oldProjResult[j];
        energy = oldProjEnergy;
    }
    free(res);
   if(flagMultiResAlsoZ == 1 || energy < oldEnergy)
        return 0;
    return -1;
}

void getNeighbourIndex(int **gridBordersIndex, int **neighboursIndex, int *noNeighbours, int n, int *gL1, int *gL2) {
    int nIndex1[2];
    int nIndex2[2];
    int gIndex1[2];
    int gIndex2[2];
    int checkCount, r1, r2;
    
    if(n == 1) {
        /* pos = [0 -1] */
        *noNeighbours = gL1[1] - gL1[0] + 1;
        nIndex1[0] = gL1[0]; nIndex1[1] = gL1[1];
        nIndex2[0] = gL2[0] - 1;
        nIndex2[1] = gL2[0] - 1;
    }
    else if(n == 2) {
        /* pos = [-1 0] */
        *noNeighbours = gL2[1] - gL2[0] + 1;
        nIndex2[0] = gL2[0]; nIndex2[1] = gL2[1];
        nIndex1[0] = gL1[0] - 1;
        nIndex1[1] = gL1[0] - 1;
    }
    else if(n == 3) {
        /* pos = [0 1] */
        *noNeighbours = gL1[1] - gL1[0] + 1;        
        nIndex1[0] = gL1[0]; nIndex1[1] = gL1[1];
        nIndex2[0] = gL2[1] + 1;
        nIndex2[1] = gL2[1] + 1;
    }    
    else if(n == 4) {
        /* pos = [1 0] */
        *noNeighbours = gL2[1] - gL2[0] + 1;        
        nIndex2[0] = gL2[0]; nIndex2[1] = gL2[1];
        nIndex1[0] = gL1[1] + 1;
        nIndex1[1] = gL1[1] + 1;        
    }
    else
        printf("getNeighbourIndex: Unknown neighbour\n");
    gIndex1[0] = nIndex1[0]; gIndex1[1] = nIndex1[1];
    gIndex2[0] = nIndex2[0]; gIndex2[1] = nIndex2[1];
    if(n == 1)      { gIndex2[0] = gIndex2[1] = nIndex2[0] + 1; }
    else if(n == 2) { gIndex1[0] = gIndex1[1] = nIndex1[0] + 1; }
    else if(n == 3) { gIndex2[0] = gIndex2[1] = nIndex2[0] - 1; }
    else if(n == 4) { gIndex1[0] = gIndex1[1] = nIndex1[0] - 1; }
    /* Check if neighbours are out of bounds */
    if(nIndex1[0] < 1 || nIndex2[0] < 1 ||
            nIndex1[1] > imageSize1 || nIndex2[1] > imageSize2) {
                *noNeighbours = 0;
                *neighboursIndex = NULL;
                *gridBordersIndex = NULL;
                return;
    }
    /* Allocate memory and return the real indices */
    /* for the grid members */
    *gridBordersIndex = (int *) malloc((*noNeighbours) * sizeof(int));
    checkCount = 0;
    for(r1 = gIndex1[0]; r1 <= gIndex1[1]; r1++)
        for(r2 = gIndex2[0]; r2 <= gIndex2[1]; r2++)
        *(*gridBordersIndex + checkCount++) = acc(r1, r2, 1);
    if(checkCount != *noNeighbours)
        printf("getNeighbourIndex: neighbours count inconsistent\n");
    /* And for neighbours */
    *neighboursIndex = (int *) malloc((*noNeighbours) * sizeof(int));
    checkCount = 0;
    for(r1 = nIndex1[0]; r1 <= nIndex1[1]; r1++)
        for(r2 = nIndex2[0]; r2 <= nIndex2[1]; r2++)
        *(*neighboursIndex + checkCount++) = acc(r1, r2, 1);
    if(checkCount != *noNeighbours)
        printf("getNeighbourIndex: neighbours count inconsistent\n");
}

double solveQuad(double a, double b, double c) {
	double dd;
	double res;

  	if(a == 0) { a = -0.00001; }
	dd = SQR(b) - 4 * a * c;
	if (dd < 0)
		printf("solveQuad: %f %f %f Oooops\n", a, b, c);
	res = (-b - sqrt(dd)) / (2*a);
//     printf("%f %f %f %f\n", a, b, c, res);    
	return res;
}

void BIDProjection(double *result, double *x)
{
    int i;
    double mes;
    double zet[MAXK];
    double Eta[MAXK];
    double VertDist[MAXK];
    double nzi[MAXK];
    double maxV; int vtx;
    int notfound;
    int c;
    double normX;
    int NOO;
    double sumTzi;
    double tzi[MAXK];
    double Tempzet[MAXK];
    double ggM; int ggExcl;
    double sumAbsX;
    
    /* Check if vector has too large variates */
    findMaximum(&ggM, &ggExcl, x);
    for(sumAbsX = 0,i = 0; i < K; i++)
        sumAbsX += fabs(x[i]);
    if(ggM > 0 && sumAbsX - ggM < 1e-5) {
        for(i = 0; i < K; i++)
            result[i] = 0;
        result[ggExcl] = 1;
        return;
    }
    /**/    
    for(mes = 0,i = 0; i < K; i++)
        mes += x[i];
    for(NOO = 0,i = 0; i < K; i++) {
        zet[i] = x[i] + (1 - mes)/K;
        if(zet[i] < 0 || zet[i] > 1)
            NOO++;
    }   
    if(NOO < 1) { /* if(NOO == 0) en fait..Blekas l'ecrit comme ca */
        for(i = 0; i < K; i++)
            result[i] = zet[i];
        return;
    }
    /*********** La piece apres "else" ***************/
    for(i = 0; i < K; i++)
        Eta[i] = 1;
    notfound = 1;
    
    for(normX = 0, i = 0; i < K; i++)
        normX += SQR(x[i]);
    for(i = 0; i < K; i++)
        VertDist[i] = 1 + normX - 2*x[i];
    c = 2;
    for(i = 0; i < K; i++)
        nzi[i] = x[i];
    while(notfound > 0) {
        findMaximum(&maxV, &vtx, VertDist);
        VertDist[vtx] = -1.0;
        Eta[vtx] = 0;
        for(sumTzi = 0, i = 0; i < K; i++) {
            tzi[i] = Eta[i]*nzi[i];
            sumTzi += tzi[i];
        }
        if(K - c + 1 <= 0) {
            printf("MEX_BIDProjection: Hey..this message shouldnt show up!\n");
            return;
        }
        for(i = 0; i < K; i++)
            Tempzet[i] = tzi[i] + (1 - sumTzi)/(K - c + 1);
        for(i = 0; i < K; i++)
            Tempzet[i] = Tempzet[i] * Eta[i];

        for(NOO = 0,i = 0; i < K; i++)
            if(Tempzet[i] < 0 || Tempzet[i] > 1)
                NOO++;
        if(NOO < 1) {
            for(i = 0; i < K; i++)
                result[i] = Tempzet[i];
            notfound = 0;
        }
        c++;
    }
    return;
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

int getNeighbourInfo(double *pos, int n)
{
    if(n == 1) {
        pos[0] = 0; pos[1] = -1;
        return 1;
    }
    if(n == 2) {
        pos[0] = -1; pos[1] = 0;
        return 2;
    }
    if(n == 3) {
        pos[0] = 0; pos[1] = 1;
        return 1;
    }
    if(n == 4) {
        pos[0] = 1; pos[1] = 0;
        return 2;
    }
}

long acc(int dir1, int dir2, int dir3) {
    /* Dir1: ranges     1..imageSize1
     * Dir2: ranges     1..imageSize2
     * Dir3: ranges     1..K
     *
     * Result ranges    0..N*K -1
     */
    long a;
    a = ((dir1-1) + (dir2-1)*imageSize1 + (dir3-1)*imageSize1*imageSize2);
//     if(a < 0 || a > imageSize1*imageSize2*K - 1)
//         printf("acc: Value %d out of bounds\n", a);
	return a;
}

long accB(int dir1, int dir2) {
    long a;
    a = ((dir1-1) + (dir2-1)*DIRECTIONS);
//     if(a < 0 || a > DIRECTIONS*K - 1)
//         printf("accB: Value %d out of bounds\n", a);    
	return a;
}

long accU(int dir1, int dir2, int dir3) {
    /* Dir1: ranges         1..Neighbour (usu.4)
     * Dir2: ranges         1..K
     * Dir3: ranges         0..N-1
     *
     * Result ranges        0..
     */
    long a;
    a = ((dir1-1) + (dir2-1)*NHOODSIZE + (dir3)*NHOODSIZE*K);
//     if(a < 0 || a > NHOODSIZE*K*imageSize1*imageSize2 - 1)
//         printf("acc: Value %d out of bounds\n", a);
	return a;
}
