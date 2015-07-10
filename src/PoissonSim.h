/**
 * Header file containing definitions of structs and function declarations used to simulate from spatial poisson model
 * 
 * Author: ggilani, 27/03/15
 */

#include "math.h"

#define MAX_PATCH 100
#define RADIUS 6371
#define PI 3.1415926535

// parameter struct
typedef struct PARAMS {
	int nPatch, tmax, nruns; //added nruns to be able to do multiple different runs at the same time
	double r0, pKernel;
	double genTimeShape, genTimeRate;
	int genTimeMax;
	double *genTime;
	//double *distMat[MAX_PATCH];
	long seed1, seed2; //seeds for initialising random numbers
} params;

// patch struct
typedef struct PATCH {
	int N, S; // number of susceptibles, total number of people in patch
	double lat, lon; //co-ordinates of patch centroid
	double infectiousness, forceInf, repRate;
	double *kernel;
	int *true_inc, *obs_inc;
} patch;

//function declarations
void initParams(int , char *[], params *, int *, int*, char *, char *);
void readInPatchInfo(char *, patch *, int);
void initGenTime(params *);
void initPatches(params *, patch *, char *);
void resetPatches(params *, patch *);
void updateInfectiousness(patch *, params *, int);
void updateForceOfInf(patch *, params *);
double expKernel(double, double);
void initEpidemic(int, int, int, patch *, params *);
void writeResults(char *, patch *, params *, int);
void updateCases(int, patch *, params*);
double getPatchDist(patch, patch);
double degreesToRad(double);
double gamm(double);
double gammaPdf(double,double,double);