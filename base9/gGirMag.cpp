#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "evolve.h"
#include "linInterp.h"
#include "binSearch.h"
#include "gGirMag.h"

#define AGELO  0
#define AGEHI  1
#define MASSLO 0
#define MASSHI 1
#define FEHLO  0
#define FEHHI  1

extern int    verbose, useFilt[FILTS], needMassNow;
extern double globalMags[FILTS];
extern double ageLimit[2];

static int  gNumEntries[N_GIR_Z], gBoundary[2][N_GIR_Z][N_GIR_AGES];
static double gIsoMass[2][N_GIR_Z][MAX_GIR_ENTRIES],
              gIsoMag[N_GIR_Z][N_GIR_FILTS][MAX_GIR_ENTRIES], 
              gLogAge[N_GIR_Z][N_GIR_AGES], gFeH[N_GIR_Z],
              gAGBt[N_GIR_Z][N_GIR_AGES];

// Set by derive_AGBt_zmass and used by getGirardiMags
static int iAge, iFeH;
static double currentAge, currentFeH;

static char tempFile[100];
static void getFileName(char *path, int z, int filterSet);


void loadGirardi(char *path, int filterSet){

  /****************************************************************************
   ***********************    filterSet = UBVRIJHK      ***********************
   *********************** Girardi models, standard mags **********************
   ** For simplicity, the 8 input files are assumed to have 50 ages each,    **
   ** which is the minimum set among all of these isochrones.  These ages    **
   ** span the range log(age) = 7.80 to 10.25.                               **
   ** The eight files and their metallicities are                            **
   **   iso_stan_z0.50.dat     Z = 0.0                                       **
   **   iso_stan_z0001.50.dat  Z = 0.0001                                    **
   **   iso_stan_z0004.50.dat  Z = 0.0004                                 **
   **   iso_stan_z001.50.dat   Z = 0.001                                  **
   **   iso_stan_z004.50.dat   Z = 0.004                                  **
   **   iso_stan_z008.50.dat   Z = 0.008                                  **
   **   iso_stan_z019.50.dat   Z = 0.019                                  **
   **   iso_stan_z030.50.dat   Z = 0.030                                  **
   ***************************************************************************/

  /****************************************************************************
   **********************        filterSet == ACS       ***********************
   ********************** Girardi models, ACS mags ****************************
   ** For simplicity, the input files are assumed to have 50 ages, just as   **
   ** the above-used Girardi standard mag isochrones.  This meant reducing   **
   ** the 50-74 ages by limiting them to log(age) = 7.80 to 10.25, which is  **
   ** fine for now as the rejected isochrones are younger than any of the    **
   ** clusters for which we have data.  The metallicities are also all the   **
   ** same.  These are the same isochrones, just convolved by Girardi with   **
   ** different filters.  The file names are different only in that the      **
   ** "stan" becomes "acs".                                                  **
   ***************************************************************************/
  
  int z, a,i,filt;
  double thisLogAge,lastLogAge;
  FILE *pGirardi;
  char line[240];

  if(filterSet != UBVRIJHK && filterSet != SDSS && filterSet != ACS){
    printf("\nFilter set %d not available on Girardi models.  Exiting...\n",filterSet);
    exit(1);
  }  
  
  for(z=0 ; z < N_GIR_Z ; z++) {                                      // foreach Girardi metallicity/isochrone file 
    
    for(a=0 ; a < N_GIR_AGES ; a++) {                                 // initialize age/boundary pointers 
      gBoundary[LOW][z][a]  = 0;
      gBoundary[HIGH][z][a] = 0;
      gLogAge[z][a]         = 0.0;
    }
    
    for(i=0;i < MAX_GIR_ENTRIES;i++) {                                // initialize array of model parameters 
      gIsoMass[ZAMS][z][i]     = 0.0;                                 // ZAMS mass
      gIsoMass[NOW][z][i] = 0.0;                                      // actual mass of star after mass loss 
      for(filt=0;filt < N_GIR_FILTS;filt++) gIsoMag[z][filt][i] = 99.999;   // U through K absolute mags 
    }
    
    getFileName(path, z, filterSet);
    //fscanf(pModelList,"%s",tempFile);                                 // work on one Girardi model at a time 
    if((pGirardi = fopen(tempFile,"r")) == NULL) {                    // open file 
      printf("\n\n file %s was not found - exiting\n",tempFile);
      exit(1);
    }
    
    
    i          = 0;
    a          = 0;                                                   // a = [0,49], gLogAge contains the 50 ages 
    lastLogAge = 0.0;
    while(fgets(line,240,pGirardi) != NULL && i < MAX_GIR_ENTRIES) {  // load first Z=# Girardi model for all ages 
      if(line[0] != '#') {
        if(filterSet == UBVRIJHK) {                                     // Girardi UBVRIJHK isocrhones 
          sscanf(line,"%lf %lf %lf %*f %*f %*f %*f %lf %lf %lf %lf %lf %lf %lf %lf %*f",
                 &thisLogAge,&gIsoMass[ZAMS][z][i],&gIsoMass[NOW][z][i],&gIsoMag[z][0][i],&gIsoMag[z][1][i],
                 &gIsoMag[z][2][i],&gIsoMag[z][3][i],&gIsoMag[z][4][i],&gIsoMag[z][5][i],&gIsoMag[z][6][i],
                 &gIsoMag[z][7][i]);
        }
        // for filterSet == ACS, use same set of variables but now have F435W F475W F550M F555W F606W F625W F775W F814W 
        // absolute mags, instead of UBVRIJHK absolute mags 
        if(filterSet == ACS) {                              // Girardi hST/ACS/WF isochrones 
          sscanf(line,"%lf %lf %lf %*f %*f %*f %*f %lf %lf %lf %lf %lf %lf %*f %*f %lf %lf %*f %*f %*f",
                 &thisLogAge,&gIsoMass[ZAMS][z][i],&gIsoMass[NOW][z][i],&gIsoMag[z][0][i],&gIsoMag[z][1][i],
                 &gIsoMag[z][2][i],&gIsoMag[z][3][i],&gIsoMag[z][4][i],&gIsoMag[z][5][i],&gIsoMag[z][6][i],
                 &gIsoMag[z][7][i]);
        }
        if(fabs(lastLogAge - thisLogAge) > EPS) {             // find model boundaries for ease/speed later 
          gBoundary[LOW][z][a]  = i;                          // gBoundary contains beginning boundary of age 
          gLogAge[z][a] = thisLogAge;
          if(a > 0){
            gAGBt[z][a-1] = gIsoMass[ZAMS][z][i-1];
            gBoundary[HIGH][z][a-1] = i-1;
          }
          
          lastLogAge        = thisLogAge;
          a++;
        }
        i++;
      }
    }
    gAGBt[z][a-1] = gIsoMass[ZAMS][z][i-1];      // Add last entry to AGBt table                                                      
    gBoundary[HIGH][z][a-1] = i-1;
    gNumEntries[z] = i;
    
    
  }
  
  gFeH[0] = -5.0;                             // Girardi isochrone metallicities: close enough for Z=0.0 
  gFeH[1] = -2.278754;                        //                                  equiv to Z=0.0001 
  gFeH[2] = -1.676694;                        //                                             0.0004 
  gFeH[3] = -1.278754;                        //                                             0.001 
  gFeH[4] = -0.676694;                        //                                             0.004 
  gFeH[5] = -0.375664;                        //                                             0.008 
  gFeH[6] = 0.0;                              // (i.e., Z_solar=0.019)                       0.019
  gFeH[7] = 0.198368;                         //                                             0.030 
  
  ageLimit[0] = gLogAge[0][0];
  ageLimit[1] = gLogAge[0][N_GIR_AGES - 1];

  fclose(pGirardi);
  
}

static void getFileName(char *path, int z, int filterSet){
	
	char fileNames[][5]={"0\0","0001\0","0004\0","001\0","004\0","008\0","019\0","030\0"};
	
	strcpy(tempFile,"\0");
	strcat(tempFile, path);
	if(filterSet == ACS) strcat(tempFile, "gIsoACS/iso_acs_z\0");
	else strcat(tempFile, "gIsoStan/iso_stan_z\0");
	strcat(tempFile, fileNames[z]);
	strcat(tempFile, ".50.dat\0");
	
}


double deriveGirAgbTip(double newFeH, double newY, double newAge)

/****************************************************************************************
last update: 12nov07

Derive AGBt mass (actually the ZAMS mass for the appropriate AGBt star) for a given 
cluster age, interpolating in isochrones as necessary.  
****************************************************************************************/

{
	
	double AGBt_zmass, AGBt_zmass_lo, AGBt_zmass_hi;
	double interpAge[2], interpFeH[2];
	
	int    binarySearch (double *searchArray, int size, double searchItem);
	
	if(newAge < 7.80) {
		if(verbose) printf("\n Requested age too young. (drv_g_AGB_m.c)");
		return 0.0;
	}
	if(newAge > 10.25) {
		if(verbose) printf("\n Requested age too old. (drv_g_AGB_m.c)");
		return 0.0;
	}
	if(newFeH < gFeH[0]){
		if(verbose) printf("\n Requested FeH too low. (drv_g_AGB_m.c)");
		return 0.0;
	}
	if(newFeH > gFeH[N_GIR_Z - 1]){
		if(verbose) printf("\n Requested FeH too high. (drv_g_AGB_m.c)");
		return 0.0;
	}
	
	interpAge[HIGH] = ceil(20. * newAge) / 20.;			// round up to nearest 0.05 
	interpAge[LOW]  = interpAge[HIGH] - 0.05;
	
	iAge            = (int)(rint((interpAge[LOW]-7.8)*20));
	
	if(verbose == 2) {
		printf("\n For AGBt: Using log(age)=%.2f bounded by interpAge[HIGH] = %.3f interpAge[LOW] = %.3f\n",
			   newAge,interpAge[HIGH],interpAge[LOW]);
	}
	
	//Find metallicity boundaries
	iFeH   = binarySearch(gFeH, N_GIR_Z, newFeH); // function returns lower bound 
	
	interpFeH[LOW]  = gFeH[iFeH];
	interpFeH[HIGH] = gFeH[iFeH + 1];
	
	AGBt_zmass_lo = linInterp(interpAge[LOW], interpAge[HIGH], gAGBt[iFeH][iAge], gAGBt[iFeH][iAge + 1], newAge);
	AGBt_zmass_hi = linInterp(interpAge[LOW], interpAge[HIGH], gAGBt[iFeH + 1][iAge], gAGBt[iFeH + 1][iAge + 1], newAge);
	
	AGBt_zmass    = linInterp(interpFeH[LOW], interpFeH[HIGH], AGBt_zmass_lo, AGBt_zmass_hi, newFeH);
	
	currentAge = newAge;
	currentFeH = newFeH;
	
	return AGBt_zmass;
	
}


double getGirardiMags(double zamsMass)

/****************************************************************************************
last update: 20jul10

Derive UBVRIJHK mags via 3-D interpolating in FeH, mass, and age within the Girardi isochrones
(for log(Age) = 7.80 to 10.25) for main sequence and red giant branch stars.  Model lists are 
in increasing mass order and all have the same low mass value.

Uses static variables iAge and iFeH.
****************************************************************************************/

{
  int    a, f, filt;
  double massNow=0.0;
  double modelAge[2], ageMassNow[2][2], ageMag[2][2][N_GIR_FILTS];
  double modelFeH[2], fehMassNow[2], fehMag[2][N_GIR_FILTS];

  modelAge[AGELO] = gLogAge[0][iAge];
  modelAge[AGEHI] = gLogAge[0][iAge + 1];
  modelFeH[FEHLO] = gFeH[iFeH];
  modelFeH[FEHHI] = gFeH[iFeH + 1];

  //logAge = 0;//currentAge;
  //FeH = 0;//currentFeH;

  // 3-D (FeH, mass. and age) interpolate U through K mags and (current, not zams) masses 
  // (rectilinear interpol requires moving along constant axes first, and it is 
  //  FeH and age that are uniform between grids. so interpolate first in mass, then age, then FeH)
  for(f=0;f<2;f++){
    // Interpolate in Mass for each Age
    for(a=0;a<2;a++) ageMassNow[f][a] = interpInMass(iAge + a, zamsMass, iFeH + f, &(ageMag[f][a][0])); 
    // Interpolate in Age for each FeH
    for(filt=0;filt < N_GIR_FILTS;filt++) {
      if(useFilt[filt] == 0) continue;			// for speed, don't bother with this calculation 
      fehMag[f][filt] = linInterp(modelAge[AGELO], modelAge[AGEHI], ageMag[f][AGELO][filt], ageMag[f][AGEHI][filt], currentAge);
    }
    if(needMassNow) {					// simCluster.c needs the resulting number 
      fehMassNow[f]      = linInterp(modelAge[AGELO], modelAge[AGEHI], ageMassNow[f][AGELO], ageMassNow[f][AGEHI], currentAge);
    }
    else fehMassNow[f] = 0.0;				// for speed, mcmc doesn't care about this value 
  }

  // Interpolate in FeH
  for(filt=0;filt < N_GIR_FILTS;filt++) {
    if(useFilt[filt] == 0) continue;			// for speed, don't bother with this calculation 
    globalMags[filt]  = linInterp(modelFeH[FEHLO],modelFeH[FEHHI],fehMag[FEHLO][filt],fehMag[FEHHI][filt],currentFeH);
  }
  if(needMassNow) {					// simCluster.c needs the resulting number 
    massNow = linInterp(modelFeH[FEHLO],modelFeH[FEHHI],fehMassNow[FEHLO],fehMassNow[FEHHI],currentFeH);
  }

  return massNow;

}

// Helper function for getGirardiMags
double interpInMass(int whichAgeIndex, double zamsMass, int whichFeHIndex, double *ageMag){

  int    filt, lo ,mid, hi, i;
  double modelMass[2][2], modelMag[2][N_GIR_FILTS];
  double tempMass = 0.0;

  lo = gBoundary[LOW][whichFeHIndex][whichAgeIndex];             // start at lowest and highest masses for this FeH and age
  hi = gBoundary[HIGH][whichFeHIndex][whichAgeIndex];
  i = -1;                                                // If this doesn't get set below, code will have a segmentation violation
                                                         // and you'll know something's wrong (this should never happen)
 
  if(zamsMass >= gIsoMass[ZAMS][whichFeHIndex][hi]) {		 // if zamsMass larger than largest value, use largest
     for(filt=0;filt < N_GIR_FILTS;filt++){                    // (this can happen when the mass is higher than the ABGt
       if(useFilt[filt]) ageMag[filt] = gIsoMag[whichFeHIndex][filt][hi]; //   for this particular age)and FeH)
     }
    if(needMassNow) tempMass = gIsoMass[ZAMS][whichFeHIndex][hi];
    return tempMass;
  }
  
  else if (zamsMass <= gIsoMass[ZAMS][whichFeHIndex][lo]) {	// if zamsMass smaller than smallest value, use smallest
    i = lo;
  }
  
  else {
    while(1) {	// binary search on zamsMass
      mid = ((lo + hi) >> 1);
      
      if(gIsoMass[ZAMS][whichFeHIndex][mid-1] <= zamsMass && zamsMass <= gIsoMass[ZAMS][whichFeHIndex][mid]) { // found upper bounding mass
        i=mid-1;
        break;
      }
      
      if(lo >= hi) {
        printf("ERROR: BINARY SEARCH FAILURE gGirmag\n");
        break;
      }
      
      if(zamsMass > gIsoMass[ZAMS][whichFeHIndex][mid]) lo = mid + 1;
      if(zamsMass < gIsoMass[ZAMS][whichFeHIndex][mid]) hi = mid - 1;
    }
  }
  
  modelMass[MASSLO][ZAMS] = gIsoMass[ZAMS][whichFeHIndex][i];
  modelMass[MASSHI][ZAMS] = gIsoMass[ZAMS][whichFeHIndex][i+1];
  modelMass[MASSLO][NOW]  = gIsoMass[NOW][whichFeHIndex][i];
  modelMass[MASSHI][NOW]  = gIsoMass[NOW][whichFeHIndex][i+1];
  for(filt=0;filt < N_GIR_FILTS;filt++) {
    if(useFilt[filt] == 0) continue;			// for speed, don't bother with this calculation 
    modelMag[MASSLO][filt]  = gIsoMag[whichFeHIndex][filt][i];
    modelMag[MASSHI][filt]  = gIsoMag[whichFeHIndex][filt][i+1];
  }
  
  for(filt=0;filt < N_GIR_FILTS;filt++) {			// grids. so move along age, i.e. interpolate first in mass) 
    if(useFilt[filt] == 0) continue;			// for speed, don't bother with this calculation 
    ageMag[filt] = linInterpExtrap(modelMass[MASSLO][ZAMS], modelMass[MASSHI][ZAMS], 
                             modelMag[MASSLO][filt], modelMag[MASSHI][filt], zamsMass);
  }
  if(needMassNow){
    tempMass = linInterpExtrap(modelMass[MASSLO][ZAMS], modelMass[MASSHI][ZAMS], 
                                 modelMass[MASSLO][NOW], modelMass[MASSHI][NOW], zamsMass);
  }

  return tempMass;

}



int reverseBinarySearch (double *searchArray, int size, double searchItem);

double wdPrecLogAgeGir(double thisFeH, double thisY, double zamsMass)

/*************************************************************************************
last update: 12nov07

Determine WD precursor age by 2-D interpolating among the AGBt mass versus age values.
Note that the appropriate AGBt mass and lifetime is not the ZAMS mass and lifetime of 
the star currently at the AGBt, but rather refers to the properties of the potentially
higher mass and younger AGBt star that was the WD precursor.
*************************************************************************************/

{


  int    thisIndexAge[2], f;
  double logAge[2], AGBt_zmass[2],wdPrecLogAge[2], FeHLo, FeHHi,temp;

  AGBt_zmass[HIGH]=AGBt_zmass[LOW]=0.0;

  FeHLo = gFeH[iFeH];
  FeHHi = gFeH[iFeH+1];

  // Find wdPrecLogAge for the lower and upper metallicity cases
  for(f=0;f<2;f++){
    
    if(zamsMass < gAGBt[iFeH + f][N_GIR_AGES-1]) {			// possible if cluster older than logAge=9.0
      wdPrecLogAge[f] = gLogAge[iFeH + f][N_GIR_AGES-1];		        // FOR NOW just use logAge = 9.0 
      if(verbose) printf(" %.3f Mo < smallest AGBt (%.2f) model mass.  Setting precursor log age to %.3f.\n",
                         zamsMass,gAGBt[iFeH + f][N_GIR_AGES-1],wdPrecLogAge[f]);
    }
    else if (zamsMass > gAGBt[iFeH + f][0]){
      wdPrecLogAge[f] = -2.7 * log10(zamsMass/gAGBt[iFeH + f][0]) + gLogAge[iFeH + f][0];
      if(verbose) printf(" %.3f Mo > largest AGBt (%.2f) model mass.  Extrapolating precursor log age.\n",
                         zamsMass,gAGBt[iFeH + f][0]);
    }

    else{
      thisIndexAge[LOW]    = reverseBinarySearch(gAGBt[iFeH + f],N_GIR_AGES,zamsMass);  // Because masses are in reverse order 
      thisIndexAge[HIGH]   = thisIndexAge[LOW] + 1;
      
      logAge[LOW]      = gLogAge[iFeH + f][thisIndexAge[LOW]];
      logAge[HIGH]     = gLogAge[iFeH + f][thisIndexAge[HIGH]];
      
      AGBt_zmass[LOW]  = gAGBt[iFeH + f][thisIndexAge[LOW]];
      AGBt_zmass[HIGH] = gAGBt[iFeH + f][thisIndexAge[HIGH]];
      
      // Linearly interpolate in mass
      wdPrecLogAge[f] = linInterp(AGBt_zmass[LOW], AGBt_zmass[HIGH], logAge[LOW], logAge[HIGH], zamsMass);
    }
  }

  // Linearly interpolate in FeH
  temp    = linInterp(FeHLo, FeHHi, wdPrecLogAge[LOW], wdPrecLogAge[HIGH], thisFeH); 

  return temp;

}

  

