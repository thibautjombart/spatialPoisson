#include "PoissonSim.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "randlib_par.h"
#include "string.h"
//#include <boost/math/distributions/gamma.hpp>
//#include <boost/random.hpp>

//typedef boost::mt19937 RNGType;
//RNGType rng;

int main(int argc, char *argv[])
{
	params P;
	patch *patches;
	int t;
	int nCases, patchID, k;
	char outFileName[1024], patchFileName[1024];


	//start by initialising parameter struct
	initParams(argc, argv, &P, &nCases, &patchID, outFileName, patchFileName);

	//space for random numbers
	initSeeds(P.seed1,P.seed2);	

	//allocate memory for number of patches
	if(!(patches=(patch *)malloc(P.nPatch*sizeof(patch))))
	{
		fprintf(stderr,"Unable to allocate enough storage for patches!\n");
		exit(-1);
	}
	//initialise patches
	initPatches(&P, patches, patchFileName);

	for(k=0;k<P.nruns;k++)
	{
		//reset time to zero
		t=0;

		//initialise the epidemic
		initEpidemic(t, nCases, patchID, patches, &P);

		//start main loop
		for(t=1;t<P.tmax;t++)
		{
			//update infectiousness
			updateInfectiousness(patches,&P,t);
			//update force of infection on each patch
			updateForceOfInf(patches,&P);
			//draw new true and observed incidences
			updateCases(t,patches,&P);
		}

		//write out results to file
		writeResults(outFileName,patches,&P,k);

		//reset patches
		resetPatches(&P, patches);

	}


}

void initParams(int argc, char *argv[], params *P, int *nCases, int *patchID, char *outFileName, char *patchFileName)
{

	int i;
	int gotNPatch,gotTMax,gotR0,gotKernelParam,gotGenTime,gotGenParams,gotPatchFile,gotInitCases,gotInitPatch,gotOutFile,gotNRuns;
	//default parameters in case some things are specified:
	int defCases=1;
	int defPatch=0;
	int defTMax=100;
	int defRuns=1; //default number of runs
	double defR0=1.1;
	//double defRepRate=1; //in default case, all cases are reported
	double defKernel=0.015; //combined with default of R0 of 1.1, gives effective R0 of ~0.5, 0.25, 0.1 for distances of 50, 100, 150km
	int defGenTime=10; //truncation of generation time distribution
	double defGenTimeShape=1.5;
	double defGenTimeRate=0.2;

	//early error warning if not enough parameters supplied - should always have nPatch, patchFile, outFile, 2 seeds
	if(argc<8)
	{
		fprintf(stderr,"Not enough parameters supplied. Correct usage: PoissonSim.exe -nPatch N -patchFile patchFileString -outFile outFileString [-tMax T -R0 R0 -kernelParam k -nruns n -genTimeMax t -genTimeParams shape rate -initialCases nCases -initialPatch P] seed1 seed2\n");
	}

	gotNPatch=gotTMax=gotR0=gotKernelParam=gotGenTime=gotGenParams=gotPatchFile=gotInitCases=gotInitPatch=gotOutFile=gotNRuns=0;
	
	//get seeds for random number initialisation
	i=argc-2;
	sscanf(argv[i],"%li",&(P->seed1));
	sscanf(argv[i+1],"%li",&(P->seed2));

	//initialise parameters from the command line
	for(i=1;i<argc-2;i++)
	{
		if(strcmp(argv[i],"-nPatch")==0)
		{
			sscanf(argv[i+1],"%i",&(P->nPatch));
			gotNPatch=1;
		}
		else if(strcmp(argv[i],"-tMax")==0)
		{
			sscanf(argv[i+1],"%i",&(P->tmax));
			gotTMax=1;
		}
		else if(strcmp(argv[i],"-R0")==0)
		{
			sscanf(argv[i+1],"%lg",&(P->r0));
			gotR0=1;
		}
		else if(strcmp(argv[i],"-nruns")==0)
		{
			sscanf(argv[i+1],"%i",&(P->nruns));
			gotNRuns=1;
		}
		else if(strcmp(argv[i],"-kernelParam")==0)
		{
			sscanf(argv[i+1],"%lg",&(P->pKernel));
			gotKernelParam=1;
		}
		else if(strcmp(argv[i],"-genTimeMax")==0)
		{
			sscanf(argv[i+1],"%i",&(P->genTimeMax));
			gotGenTime=1;
		}
		else if(strcmp(argv[i],"-genTimeParams")==0)
		{
			sscanf(argv[i+1],"%lg",&(P->genTimeShape));
			sscanf(argv[i+2],"%lg",&(P->genTimeRate));
			gotGenParams=1;
		}
		else if(strcmp(argv[i],"-patchFile")==0)
		{
			sscanf(argv[i+1],"%s",patchFileName);
			gotPatchFile=1;
		}
		else if(strcmp(argv[i],"-initialCases")==0)
		{
			sscanf(argv[i+1],"%i",nCases);
			gotInitCases=1;
		}
		else if(strcmp(argv[i],"-initialPatch")==0)
		{
			sscanf(argv[i+1],"%i",patchID);
			gotInitPatch=1;
		}
		else if(strcmp(argv[i],"-outFile")==0)
		{
			sscanf(argv[i+1],"%s",outFileName);
			gotOutFile=1;
		}

	}
	//set default values if needed
	if(gotTMax==0)
	{
		fprintf(stderr,"Max simulation time has not been specified. Using default value of tmax: %i\n",defTMax);
		P->tmax=defTMax;
	}
	if(gotR0==0)
	{
		fprintf(stderr,"R0 has not been specified. Using default value of R0: %lg\n",defR0);
		P->r0=defR0;
	}
	if(gotNRuns==0)
	{
		fprintf(stderr,"No. of runs has not been specified. Using default value of nruns: %i\n",defRuns);
		P->nruns=defRuns;
	}
	if(gotKernelParam==0)
	{
		fprintf(stderr,"Kernel parameter has not been specified. Using default value of: %lg\n",defKernel);
		P->pKernel=defKernel;
	}
	if(gotGenTime==0)
	{
		fprintf(stderr,"Max generation time has not been specified. Using default value of: %i\n",defGenTime);
		P->genTimeMax=defGenTime;
	}
	if(gotGenParams==0)
	{
		fprintf(stderr,"Generation time parameters have not been specified. Using default value of shape: %lg, scale: %lg\n",defGenTimeShape,defGenTimeRate);
		P->genTimeShape=defGenTimeShape;
		P->genTimeRate=defGenTimeRate;
	}
	if(gotInitCases==0)
	{
		fprintf(stderr,"Initial number of cases has not been specified. Using default value of: %i\n",defCases);
		*nCases=defCases;
	}
	if(gotInitPatch==0)
	{
		fprintf(stderr,"Initialisation patch has not been specified. Using default value of: %i\n",defPatch+1);
		*patchID=defPatch;
	}
	//notify errors to do with number of patches and distance matrix file
	if(gotNPatch==0)
	{
		fprintf(stderr,"Number of patches hasn't been specified - it must match the number given in the distance matrix file\n");
		exit(-2);
	}
	if((P->nPatch>1)&&(gotPatchFile==0))
	{
		fprintf(stderr,"File containing distance matrix hasn't been specified\n");
		exit(-3);
	}
	if(gotOutFile==0)
	{
		fprintf(stderr,"Output file name hasn't been specified\n");
		exit(-4);
	}

	//initialise generation time vector
	initGenTime(P);
}

void initGenTime(params *P)
{
	int i;
	double genTimeSum=0.0;
	//boost::math::gamma_distribution<> gammadist(P->genTimeShape,P->genTimeScale);

	//allocate memory for generation time vector
	if(!(P->genTime=(double *)malloc(P->genTimeMax*sizeof(double))))
	{
		fprintf(stderr,"Unable to allocate enough storage for generation time vector!\n");
		exit(-5);
	}

	//calculate generation time distribution
	for(i=0;i<P->genTimeMax;i++)
	{
		P->genTime[i]=gammaPdf((double)i+1,P->genTimeShape,P->genTimeRate); //index zero stores day 1, index genTimeMax-1 stores day genTimeMax
		genTimeSum+=P->genTime[i];
	}

	//normalise generation time distribution
	for(i=0;i<P->genTimeMax;i++)
	{
		P->genTime[i]/=genTimeSum;
	}
}

void initPatches(params *P, patch *patches, char *patchFileName)
{
	int i,j;
	double kernelsum,dist;

	//allocate memory for each patch to store true incidence, obs incidence and kernel
	for(i=0;i<P->nPatch;i++)
	{
		if(!(patches[i].kernel=(double *) malloc(P->nPatch*sizeof(double))))
		{
			fprintf(stderr,"Unable to allocate enough storage for patch kernel vector!\n");
			exit(-6);
		}
		if(!(patches[i].true_inc=(int *) malloc(P->tmax*sizeof(int))))
		{
			fprintf(stderr,"Unable to allocate enough storage for true incidence!\n");
			exit(-7);
		}
		if(!(patches[i].obs_inc=(int *) malloc(P->tmax*sizeof(int))))
		{
			fprintf(stderr,"Unable to allocate enough storage for observed incidence!\n");
			exit(-8);
		}
	}

	//read in 
	readInPatchInfo(patchFileName, patches, P->nPatch);

	//initialise the kernel
	for(i=0;i<P->nPatch;i++)
	{
		kernelsum=0;
		for(j=0;j<P->nPatch;j++)
		{
			dist=getPatchDist(patches[i],patches[j]);
			patches[i].kernel[j]=expKernel(dist,P->pKernel);
			kernelsum+=patches[i].kernel[j];
		}
		for(j=0;j<P->nPatch;j++)
		{
			patches[i].kernel[j]/=kernelsum;

		}
	}

	//for each patch, initialise true and observed incidence to zero, infectiousness and force of infection to zero
	for(i=0;i<P->nPatch;i++)
	{
		patches[i].forceInf=0.0;
		patches[i].infectiousness=0.0;
		for(j=0;j<P->tmax;j++)
		{
			patches[i].true_inc[j]=0;
			patches[i].obs_inc[j]=0;
		}
	}
	
}

void readInPatchInfo(char *patchFileName, patch *patches, int nPatch)
{
	int i,j,n;
	FILE *dat;
	
	//open file containing patch information
	if(!(dat=fopen(patchFileName,"r"))) 
	{
		fprintf(stderr,"Unable to open patch file\n");
		exit(-9);
	}
	//read in first value from file, which should be an int containing the number of patches included in the distance matrix file
	fscanf(dat,"%i",&n);
	//if n doesn't match the number of patches in the parameter file, there's a problem!
	if(n!=nPatch)
	{
		fprintf(stderr,"Number of patches in distance matrix file doesn't match the number of patches in the parameter struct\n");
		exit(-10);
	}

	// read in rest of patch info, in the form patch_x, patch_y, patch_N
	for(i=0;i<nPatch;i++)
	{
		fscanf(dat,"%lg %lg %i %lg",&patches[i].lat,&patches[i].lon,&patches[i].N,&patches[i].repRate);
		// set number of susceptibles equal to number of people in the patch.
		patches[i].S=patches[i].N;
	}

	//close distance matrix file
	fclose(dat);
}

double getPatchDist(patch patch1, patch patch2)
{
	double lat1,lat2,lon1,lon2,lat_diff,lon_diff,cent_angle;

	//assign lats and lons
	lat1=degreesToRad(patch1.lat);
	lon1=degreesToRad(patch1.lon);
	lat2=degreesToRad(patch2.lat);
	lon2=degreesToRad(patch2.lon);

	lat_diff=lat1-lat2;
	lon_diff=lon1-lon2;

	//calculate central angle between the two points
	cent_angle=sqrt(sin(0.5*lat_diff)*sin(0.5*lat_diff)+cos(lat1)*cos(lat2)*sin(0.5*lon_diff)*sin(0.5*lon_diff));
	cent_angle=2*asin(cent_angle);

	return RADIUS*cent_angle;
}

double degreesToRad(double degree)
{
	return (degree/180.0)*PI;
}

void resetPatches(params *P, patch *patches)
{
	int i,j;
	
	//for each patch, initialise true and observed incidence to zero, infectiousness and force of infection to zero
	for(i=0;i<P->nPatch;i++)
	{
		patches[i].forceInf=0.0;
		patches[i].infectiousness=0.0;
		// reset number of susceptibles in patch equal to population of the patch
		patches[i].S=patches[i].N;
		for(j=0;j<P->tmax;j++)
		{
			patches[i].true_inc[j]=0;
			patches[i].obs_inc[j]=0;
		}
	}
}

void updateInfectiousness(patch *patches, params *P, int tCurrent)
{
	int i,t,j;

	//loop over all patches
	for(i=0;i<P->nPatch;i++)
	{
		//reset infectiousness to zero for the patch
		patches[i].infectiousness=0.0;
		//loop over times in generation time distribution
		for(j=0;j<P->genTimeMax;j++)
		{
			t=tCurrent-(j+1);
			if(t>=0)
			{
				patches[i].infectiousness+=(patches[i].true_inc[t]*P->genTime[j]);
			}
		}
	}
	
}

void updateForceOfInf(patch *patches, params *P)
{
	int i,j;
	
	for(i=0;i<P->nPatch;i++)
	{
		//reset force of infection to zero
		patches[i].forceInf=0;
		for(j=0;j<P->nPatch;j++)
		{
			patches[i].forceInf+=(patches[j].infectiousness*patches[i].kernel[j]);
		}
		//finally multiply by R0
		patches[i].forceInf*=P->r0;
		//also scale by percentage of susceptibles left in patch
		patches[i].forceInf*=((double)patches[i].S/(double)patches[i].N);
	}
}

void initEpidemic(int t, int nCases, int patchID, patch *patches, params *P)
{
	// set number of cases in initial patch
	patches[patchID].true_inc[t]=nCases;
	// reduce number of susceptibles in the patch accordingly
	patches[patchID].S-=nCases;
	//draw appropriate observed incidence from binomial distribution
	//boost::random::binomial_distribution<> bin(patches[patchID].true_inc[t],P->repRate);
	//boost::variate_generator<RNGType&, boost::random::binomial_distribution<>> bin_sampler(rng, bin);
	//update observed incidence
	//patches[patchID].obs_inc[t]=bin_sampler();
	patches[patchID].obs_inc[t]=ignbin((long) patches[patchID].true_inc[t],patches[patchID].repRate);
}

double expKernel(double distance, double kernelParam)
{
	return exp(-kernelParam*distance);
}

void updateCases(int t, patch *patches, params *P)
{
	int i;

	for(i=0;i<P->nPatch;i++)
	{
		//set up random number distribution for poisson sampling
		if(patches[i].forceInf>0)
		{
			//boost::random::poisson_distribution<> poiss(patches[i].forceInf); //scale force of infection by number of susceptibles left in patch.
			//boost::variate_generator<RNGType&, boost::random::poisson_distribution<>> poiss_sampler(rng, poiss);
			//update true incidence
			//patches[i].true_inc[t]=poiss_sampler();
			patches[i].true_inc[t]=ignpoi(patches[i].forceInf);
		}
		else{
			patches[i].true_inc[t]=0;
		}

		//update number of susceptibles in patch
		if(patches[i].true_inc[t]<=patches[i].S)
		{
			patches[i].S-=patches[i].true_inc[t];
		}
		else
		{
			patches[i].true_inc[t]=patches[i].S;
			patches[i].S-=patches[i].true_inc[t];
		}

		//set up random number distribution for binomial sampling
		//boost::random::binomial_distribution<> bin(patches[i].true_inc[t],P->repRate);
		//boost::variate_generator<RNGType&, boost::random::binomial_distribution<>> bin_sampler(rng, bin);
		//update observed incidence
		//patches[i].obs_inc[t]=bin_sampler();
		patches[i].obs_inc[t]=ignbin((long)patches[i].true_inc[t],patches[i].repRate);
	}
}

void writeResults(char *outFileName, patch *patches, params *P, int k)
{
	int i,j;
	char OutFileNameTemp[1024];
	FILE *dat;

	sprintf(OutFileNameTemp,"%s_%i.csv",outFileName,k);
	
	if(!(dat=fopen(OutFileNameTemp,"w"))) 
	{
		fprintf(stderr,"Unable to open output file\n");
		exit(-11);
	}

	//read out parameters first: R0; reporting rate, kernel parameter
	fprintf(dat,"%i\n",P->nPatch);
	fprintf(dat,"%lg\n",P->r0);
	for(i=0;i<P->nPatch;i++)
	{
		fprintf(dat,"%lg, ",patches[i].repRate);
	}
	fprintf(dat,"\n");
	fprintf(dat,"%lg\n",P->pKernel);

	//read out true incidence for each patch: one row per patch
	for(i=0;i<P->nPatch;i++)
	{
		for(j=0;j<P->tmax;j++)
		{
			fprintf(dat,"%i, ",patches[i].true_inc[j]);
		}
		fprintf(dat,"\n");
	}
	//read out observed incidence for each patch: one row per patch
	for(i=0;i<P->nPatch;i++)
	{
		for(j=0;j<P->tmax;j++)
		{
			fprintf(dat,"%i, ",patches[i].obs_inc[j]);
		}
		fprintf(dat,"\n");
	}

	fclose(dat);
}

//
// Lanczos approximation to the gamma function. 
// 
// found on http://www.rskey.org/gamma.htm   
//
double gamm(double x) 
{
    double ret = (1.000000000190015 + 
                 76.18009172947146 / (x + 1) +  
                 -86.50532032941677 / (x + 2) + 
                 24.01409824083091 / (x + 3) +  
                 -1.231739572450155 / (x + 4) + 
                 1.208650973866179e-3 / (x + 5) + 
                 -5.395239384953e-6 / (x + 6));
    
    return ret * sqrt(2*PI)/x * pow(x + 5.5, x+.5) * exp(-x-5.5);
}

// the gamma distribution PDF
double gammaPdf(double x, double a, double b)
{
    if (x <= 0 || a <= 0 || b <= 0)
        return 0.0;
    else
        return exp(-x * b) * pow(x, a - 1.0) * pow(b, a) / gamm(a);
}