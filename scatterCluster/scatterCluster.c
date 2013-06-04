/*
  25aug10      - changed to output sigma instead of variance
  20jul10      - updated to work with different filter sets
  26jun08      - updated to new noise model based on scaling KPNO 4m
  28may08      - use apparent instead of absolute mag's
  28sep07      - cosmetic update
  13sep06      - update to add simulated field stars
  13jun06      - update to handle MS+MS binaries
  11dec05      - change output probability from FS prob to CS prob, i.e. 0 -> 1
  03dec05      - add filters used to header, making assumption it is just BVI for now
  12nov05      - output flag (=0) to indicate probability that object is a field star
  07oct05      - output stage of MS/RGB, WD star and N_WD for ease of compiling stats later
  28sep05      - change gasdev call to gennorm call
  28feb05      - temporary kludge - reject stars below 0.25 Mo - no age info and don't want
  to create too many objects that could hop below the 0.15 Mo model limit
  14sep04      - modify to work on binaries, not object 1; pass in brightLimit cut-off
  29feb04      - make variance realistic and magnitude dependent as per HST/ACS
  09dec03      - update gasdev call
  04nov03      - modify output format
  scatterCluster.c  30oct03  tvh - keep a subset of the output of simCluster.c and scatter the photometry.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "mt19937ar.h"
#include "evolve.h"
#include "structures.h"
#include "Settings.hpp"

unsigned long mt[NN];
unsigned long seed=0;
int           mti=NN+1;

static double mass1, mass2, phot[FILTS], exptime[FILTS], sigma[FILTS];
static int stage1, stage2, starID;

double signalToNoise(double mag, double exptime, int filt);
int readLine(FILE* filePtr);
int magCutoff(int firstFilt, double brightLimit, double faintLimit);
int stageCutoff(int isFS);
int scatterPhot(double limitSigToNoise);
double gen_norm(double mean, double std_dev);  // Returns a normal rv
int outputScatter(FILE* w_ptr, int isFS, double clusterMemberPrior);



int main(int argc, char *argv[])

{

    int    count, nr, nStars, wdCount, i, filt,
        filterSet, firstFilt, nFieldStars, isFS, isBD;
    double limitSigToNoise, brightLimit, faintLimit, clusterMemberPrior;
    char   filename[100], line[1000], aFilterName[10];
    FILE   *r_ptr, *w_ptr;

    struct Settings *settings = malloc(sizeof(struct Settings));
    settingsFromCLI(argc, argv, settings);
    if (settings->files.config)
    {
        makeSettings(settings->files.config, settings);
    }
    else
    {
        makeSettings("base9.yaml", settings);
    }

    settingsFromCLI(argc, argv, settings);

    /* printf("\n Enter simulated cluster file name : "); */
    /* scanf("%s",filename); */
    strcpy(filename, settings->files.output);
    strcat(filename, ".sim.out");
    if((r_ptr = fopen(filename,"r")) == NULL) {
        printf("\n\n file %s was not found - exiting ",filename);
        exit(1);
    }

    //Scan header line to figure out which photometry set is being used
    //for(i=0;i<29;i++) fscanf(r_ptr,"%*s ");
    //fscanf(r_ptr,"%s",aFilterName);

    for(i=0;i<2;i++) fscanf(r_ptr,"%*s ");
    fscanf(r_ptr,"%s",aFilterName);
    i=0;
    while(aFilterName[i]!='1') i++;
    aFilterName[i] = '\0';


    for(filterSet=0;filterSet<3;filterSet++){
        setFilterNames(filterSet);
        printf("%d %s %s\n",filterSet, aFilterName, getFilterName(0));
        if(strcmp(aFilterName,getFilterName(0)) == 0) break;
    }
    printf("filterSet = %d\n",filterSet);

    fgets(line,1000,r_ptr);		// remove rest of header line

    /* printf("\n Enter hours of exposure for noise model for each of "); */
    /* for(filt=0;filt<FILTS;filt++) */
    /*     printf("%s ",getFilterName(filt)); */
    /* //else            printf("\n Enter hours of exposure for noise model for each of band1 band2 ... band8"); */
    /* printf("\n                       e.g., 2.3 1.0 1.2 0. 0. 0. 0. 0."); */
    /* printf("\n                       where 0. exposure time means unused band. "); */
    /* scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", */
    /*       &exptime[0],&exptime[1],&exptime[2],&exptime[3], &exptime[4],&exptime[5],&exptime[6],&exptime[7], */
    /*       &exptime[8],&exptime[9],&exptime[10],&exptime[11],&exptime[12],&exptime[13]); */

    memcpy(exptime, settings->scatterCluster.exposures, 14 * sizeof(double));

    /* printf("\n Enter number of stars to keep, bright and faint and cut-off mags, and their filter: "); */
    /* scanf("%d %lf %lf %d",&nStars,&brightLimit, &faintLimit, &firstFilt);		// brightLimit primarily used to cut off RGB */

    nStars = settings->simCluster.nStars;
    brightLimit = settings->scatterCluster.brightLimit;
    faintLimit = settings->scatterCluster.faintLimit;
    firstFilt = settings->scatterCluster.relevantFilt;

    /* printf("\n Enter limiting signal-to-noise (e.g. 15) : "); */
    /* scanf("%lf",&limitSigToNoise); */
    limitSigToNoise = settings->scatterCluster.limitS2N;

    /* printf("\n Enter number of field stars to include (e.g. 0): "); */
    /* scanf("%d",&nFieldStars); */
    nFieldStars = settings->simCluster.nFieldStars;
    if(nFieldStars < 0) nFieldStars = 0;

    /* printf("\n Enter an integer seed: "); */
    /* scanf("%ld",&seed); */
    seed = settings->seed;

    /* printf("\n Enter output file name : "); */
    /* scanf("%s",filename); */
    strcpy(filename, settings->files.output);
    strcat(filename, ".sim.scatter");
    if((w_ptr = fopen(filename,"w")) == NULL) {
        printf("\n\n file %s not available for writing - exiting ",filename);
        exit(1);
    }
    /* printf("\n"); */



    // clusterMemberPrior = (double) nStars / (nStars + nFieldStars);
    // if(clusterMemberPrior > 0.99) clusterMemberPrior = 0.99;
    clusterMemberPrior = 1.0; // NS TEMPORARY!!!

    init_genrand(seed);

    //Output headers
    for(filt=0;filt<FILTS;filt++) fprintf(w_ptr, "%s ", getFilterName(filt));
/*    fprintf(w_ptr,"\n");
    for(filt=0;filt<FILTS;filt++){
        if(exptime[filt] > EPS) fprintf(w_ptr,"1 ");
        else                    fprintf(w_ptr,"0 ");
    }
*/
    fprintf(w_ptr,"id ");
    for(filt=0;filt<FILTS;filt++) if(exptime[filt] > EPS) fprintf(w_ptr,"%6s ",   getFilterName(filt));
    for(filt=0;filt<FILTS;filt++) if(exptime[filt] > EPS) fprintf(w_ptr,"sig%-5s ",getFilterName(filt));

    fprintf(w_ptr,"mass1 massRatio stage1 CMprior useDuringBurnIn\n");

    count   = 0;
    wdCount = 0;

    isFS    = 0;
    isBD    = 0;

    //Fix this line to read in however many bands there are
    while((nr = readLine(r_ptr)) != EOF) {
        // simCluster creates brown dwarfs with #'s starting at 10,001
        if(starID > 10000 && !isBD) {
            if(count < nStars) {
                printf(" Warning - not enough (%d) non-RGB stars in input file to keep desired number (%d) of stars.\n",count,nStars);
            }
            count = 0;
            nStars = 100000;
            isBD  = 1;
        }

        // simCluster creates field stars with #'s starting at 20,001
        if(starID > 20000 && !isFS) {
            count = 0;
            isFS  = 1;
        }

        if(magCutoff(firstFilt,brightLimit,faintLimit) == 0) continue;
        if(stageCutoff(isFS) == 0) continue;
        if(scatterPhot(limitSigToNoise) == 0) continue;
        if(mass1 < 0.25 && mass2 < 0.25  && stage1 != BD) continue;		// limit on modelSet=4 is 0.4 Mo - what to do?
        //if(stage1 == BD) continue;
        wdCount += outputScatter(w_ptr, isFS, clusterMemberPrior);
        count++;

        // If we have enough stars of this type
        if(count == nStars) {
            // If these are field stars, we're done
            if(isFS) break;
            // If not...
            else {
                // Skip lines in the file until you get to...
                while((nr = fscanf(r_ptr,"%d ",&starID)) != EOF) {
                    //...the brown dwarfs (which start at 10001)...
                    if(starID==10001){
                        isBD = 1;
                        break;
                    }
                    //...or the field stars (which start at 20001)
                    if(i==20001){
                        isFS = 1;
                        isBD = 1;
                        break;
                    }
                    fgets(line,1000,r_ptr);
                }
                fseek(r_ptr, -7, SEEK_CUR);

                // If this is the first field star, reset the nStars to
                // nFieldStars and read in nFieldStars in the main loop
                if(isFS){
                    nStars = nFieldStars;
                    count = 0;
                    continue;
                }
                else if(isBD){
                    nStars = 100000000;
                    count = 0;
                    continue;
                }
            }
        }
    }
    if(count < nFieldStars) {
        printf(" Warning - not enough (%d) field stars in input file to keep desired number (%d) of stars.\n",
               count,nFieldStars);
    }

    printf(" There %s %d WD%s in this scatter file, %s\n",wdCount == 1 ? "is" : "are",wdCount,wdCount ==1 ? "" : "s",filename);

    fclose(r_ptr);
    fclose(w_ptr);

    return(0);
}


static double s2nCoeffs[][2]=
{
    {9.33989,  0.3375778},	// U
    {10.0478,  0.3462758},	// B
    {10.48098, 0.368201 },	// V
    {10.71151, 0.3837847},	// R
    {10.61035, 0.3930941},	// I
    {9.282385, 0.386258 },	// J
    {9.197463, 0.3970419},	// H
    {9.024068, 0.3985604},  // K
    {9.024068, 0.3985604},  // IRAC Blue
    {9.024068, 0.3985604},  // IRAC Red
    {9.024068, 0.3985604},  // Band 1
    {9.024068, 0.3985604},  // Band 2
    {9.024068, 0.3985604},  // Band 3
    {9.024068, 0.3985604}   // Band 4
};

double signalToNoise(double mag, double exptime, int filter)

/*
  This is an approximation to the results one would obtain in one hour with the KPNO
  4m + Mosaic (UBVRI) or Flamingos (JHK) per band, assuming dark time, seeing=1.1",
  airmass=1.2, and then scaling from their by sqrt(exptime).  I further approximated
  the CCDTIME results (run on the NOAO webste) with linear fits of mag vs. log(S/N).
*/

{

    double s2n, logS2N;

    if(filter >= FILTS || filter < 0){
        printf("filter (%d) out of range - exiting\n",filter);
        exit(1);
    }

    // Scatter BD photometry at 5%
    if(filter > 7) return 1.0/0.05;

    logS2N  = s2nCoeffs[filter][0] - s2nCoeffs[filter][1] * mag;
    s2n  = pow(10., logS2N);
    s2n *= sqrt(exptime);

    return(s2n);

}

int readLine(FILE* filePtr){

    int nr,filt;
    fscanf(filePtr, "%d %lf ", &starID, &mass1);
    for(filt=0;filt<FILTS;filt++)  fscanf(filePtr, "%*f ");
    fscanf(filePtr, "%d %*f %*d %*f %*f %lf ",&stage1,&mass2);
    for(filt=0;filt<FILTS;filt++)  fscanf(filePtr, "%*f ");
    nr = fscanf(filePtr, "%d %*f %*d %*f %*f ",&stage2);
    for(filt=0;filt<FILTS;filt++)  nr = fscanf(filePtr, "%lf ",&phot[filt]);

    if(starID==EOF){
        printf("\nFatal error in readLine.  Exiting.\n");
        exit(1);
    }

    if(nr == EOF) return nr;
    return starID;

}

int magCutoff(int firstFilt, double brightLimit, double faintLimit){
    if(phot[firstFilt] > 99. && stage1 != BD) return 0;         // check if real object.  if not, skip
    if(phot[firstFilt] < brightLimit && stage1 != BD) return 0;	// typically used to remove red giants
    if(phot[firstFilt] > faintLimit && stage1 != BD) return 0;	// typically used to remove WDs and lower MS
    return 1;
}

int stageCutoff(int isFS){
    if(stage1 == NSBH || stage2 == NSBH) return 0;
    if(stage1 == WD && mass2 > 0.0 && !isFS) return 0;	// TEMPORARY KLUDGE -- ignore binaries of MS/RG + WDs and WD + WD
    if(stage2 == WD && mass1 > 0.0 && !isFS) return 0;
    return 1;
}

int scatterPhot(double limitSigToNoise){
    int filt;
    double sigToNoise;

    for(filt=(stage1==BD?8:0);filt<(stage1==BD?FILTS:8);filt++){
        if(exptime[filt] < EPS) continue;
        sigToNoise  = signalToNoise(phot[filt], exptime[filt], filt);
        if(sigToNoise < limitSigToNoise) {		// large photometric errors can lock mcmc during burnin
            printf("Warning: star %4d, mass1=%.3f, stage1=%d, filter=%d, S/N = %.3f < %.1f (user limit) - skipping.\n",
                   starID, mass1,stage1,filt,sigToNoise,limitSigToNoise);
            return 0;
        }

        sigma[filt] = 1./(sigToNoise);
        // if(sigma[filt] < 0.005) sigma[filt] = 0.005; // NS TEMPORARY!!
        if(sigma[filt] < 0.01) sigma[filt] = 0.01;
        // phot[filt] += gen_norm(0., sigma[filt]);
        phot[filt] += gen_norm(0., 3*sigma[filt]); // NS TEMPORARY!!!
    }
    return 1;
}

int outputScatter(FILE* w_ptr, int isFS, double clusterMemberPrior){
    int tempStage, filt;
    double tempMass, tempMassRatio;

    if(mass1 > mass2) {            // use to set starter mass for mcmc
        tempMass      = mass1;       // (since higher mass star dominates the photometry)
        tempMassRatio = mass2/mass1;
        tempStage     = stage1;
    }
    else {
        tempMass      = mass2;
        tempMassRatio = mass1/mass2;
        tempStage     = stage2;
    }

    //Fix outputs
    fprintf(w_ptr,"%6d ",starID);
    for(filt=0;filt<FILTS;filt++) if(exptime[filt] > EPS) fprintf(w_ptr,"%6.3f ",phot[filt]);
    for(filt=0;filt<8;filt++){
        if(exptime[filt] > EPS){
            if(stage1 == BD) fprintf(w_ptr,"%8.5f ", -1.0);
            // else fprintf(w_ptr,"%8.6f ",sigma[filt]);
            else fprintf(w_ptr,"%8.6f ",sigma[filt]*sigma[filt]); // TEMPORARY!!!
        }
    }
    for(filt=8;filt<FILTS;filt++){
        if(exptime[filt] > EPS){
            // if(stage1 == BD) fprintf(w_ptr,"%8.6f ",sigma[filt]);
            if(stage1 == BD) fprintf(w_ptr,"%8.6f ",sigma[filt]*sigma[filt]); // TEMPORARY!!!
            else fprintf(w_ptr,"%8.5f ",-1.0);
        }
    }

    fprintf(w_ptr,"%8.3f %6.3f %3d %6.3f   %d\n",tempMass,(tempMassRatio>0.001?tempMassRatio:0.00),tempStage,clusterMemberPrior,!isFS);

    if(tempStage == WD && !isFS) return 1;
    else return 0;
}



