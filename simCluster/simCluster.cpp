// last update, 26jun08

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include "mt19937ar.hpp"
#include "evolve.hpp"
#include "structures.hpp"
#include "loadModels.hpp"
#include "Settings.hpp"

double gen_norm (double mean, double std_dev);

// Used by other methods in simCluster
static int nMSRG = 0, nWD = 0, nNSBH = 0;       //, nDa=0, nDb=0;
static double wdMassTotal = 0.0, MSRGMassTotal = 0.0;

// Used by evolve.c
double ltau[2];
int aFilt = 0;

// Used by a bunch of different functions.
int verbose, needMassNow = 1, useFilt[FILTS];

// For random # generator (mt19937ar.c)
unsigned long mt[NN], seed = 0;
int mti = NN + 1;

int main (int argc, char *argv[])
{
    int i, filt, nStars, cmpnt, nFieldStars, nBrownDwarfs;
    double fractionBinary, tempU, massTotal, fractionDB, tempMod, minV, maxV, minMass = 0.15;
    char w_file[100];
    FILE *w_ptr;
    struct cluster theCluster;
    struct star theStar;

    double drawFromIMF (void);
    double genrand_res53 (void);
    void updateCount (struct star *pStar, int cmpnt);

    struct Settings *settings = malloc (sizeof (struct Settings));

    zeroSettingPointers (settings);
    settingsFromCLI (argc, argv, settings);
    if (settings->files.config)
    {
        makeSettings (settings->files.config, settings);
    }
    else
    {
        makeSettings ("base9.yaml", settings);
    }

    settingsFromCLI (argc, argv, settings);

    //macros defined in structures.h
    initCluster (&theCluster);
    initStar (&theStar);

    for (filt = 0; filt < 8; filt++)
        useFilt[filt] = 1;              // calculate all of U-K
    for (filt = 8; filt < FILTS; filt++)
        useFilt[filt] = 1;              // but not the other crap

    printf ("\n ***You are running simCluster version %.1f.***\n", VERSION);
    // clusY needed, but ignored unless modelSet = 3
    /* printf("\n Enter nSystems, WDMassUp, percentBinary, percentDB, (m-M)v, Av, logClusterAge, [Fe/H], Y, nFieldStars, nBrownDwarfs : "); */
    /* printf("(e.g., 1000 6.0 50 25 12.0 0.5 9.0 -0.3 0.27 100) : "); */
    /* scanf("%d %lf %lf %lf %lf %lf %lf %lf %lf %d %d", &theCluster.nStars, &theCluster.M_wd_up,&fractionBinary,&fractionDB, */
    /*       &theCluster.parameter[MOD],&theCluster.parameter[ABS],&theCluster.parameter[AGE], */
    /*       &theCluster.parameter[FEH],&theCluster.parameter[YYY], &nFieldStars, &nBrownDwarfs); */

    theCluster.nStars = settings->simCluster.nStars;
    theCluster.M_wd_up = settings->whiteDwarf.M_wd_up;
    fractionBinary = settings->simCluster.percentBinary;
    fractionDB = settings->simCluster.percentDB;
    theCluster.parameter[MOD] = settings->cluster.distMod;
    theCluster.parameter[ABS] = settings->cluster.Av;
    theCluster.parameter[AGE] = settings->cluster.logClusAge;
    theCluster.parameter[FEH] = settings->cluster.Fe_H;
    theCluster.parameter[YYY] = settings->cluster.Y;
    nFieldStars = settings->simCluster.nFieldStars;
    nBrownDwarfs = settings->simCluster.nBrownDwarfs;

    fractionBinary /= 100.;     // input as percentages, use as fractions
    fractionDB /= 100.;

    /* printf("\n Enter an integer seed: "); */
    /* scanf("%ld",&seed); */
    seed = settings->seed;

    /* printf("\n Run in verbose mode (0=no, 1=yes, 2=YES) ?"); */
    /* scanf("%d",&verbose); */
    verbose = settings->verbose;
    if (verbose < 0 || verbose > 2)
        verbose = 1;            // give standard feedback if incorrectly specified

    loadModels (nFieldStars, &theCluster, settings);

    if (theCluster.evoModels.mainSequenceEvol == YALE)
        minMass = 0.4;
    if (theCluster.evoModels.mainSequenceEvol == DSED)
        minMass = 0.25;

    /* printf("\n Enter CM diag output file name : "); */
    /* scanf("%s",w_file); */
    strcpy (w_file, settings->files.output);
    strcat (w_file, ".sim.out");
    if ((w_ptr = fopen (w_file, "w")) == NULL)
    {
        printf ("\n\n file %s was not available for writing - exiting ", w_file);
        exit (1);
    }
    /* printf("\n\n"); */

    nStars = 0;
    massTotal = 0.0;

    init_genrand (seed);

    //Output headers
    fprintf (w_ptr, "id  mass1 ");
    for (filt = 0; filt < FILTS; filt++)
        if (useFilt[filt])
            fprintf (w_ptr, "%s1 ", getFilterName (filt));
    fprintf (w_ptr, "stage1 wdM1 wdType1 wdLogTeff1 ltau1 mass2 ");
    for (filt = 0; filt < FILTS; filt++)
        if (useFilt[filt])
            fprintf (w_ptr, "%s2 ", getFilterName (filt));
    fprintf (w_ptr, "stage2 wdM2 wdType2 wdLogTeff2 ltau2 ");
    for (filt = 0; filt < FILTS; filt++)
        if (useFilt[filt])
            fprintf (w_ptr, "%s ", getFilterName (filt));
    fprintf (w_ptr, "\n");

    minV = 1000.0;
    maxV = -1000.0;

    ////////////////////////////////
    ///// Create cluster stars /////
    ////////////////////////////////
    // derive masses, mags, and summary stats
    for (i = 0; i < theCluster.nStars; i++)
    {                           // for all systems in the cluster
        do
        {
            theStar.U = drawFromIMF (); // create single stars in arbitrary mass order
        } while (theStar.U < minMass);
        nStars++;                       // keep track of the number of stars created
        massTotal += theStar.U;
        theStar.massRatio = 0.0;

        evolve (&theCluster, &theStar, 0);      // given inputs, derive mags for first component

        fprintf (w_ptr, "%4d %7.3f ", i + 1, getMass1 (&theStar, &theCluster)); // output primary star data
        for (filt = 0; filt < FILTS; filt++)
        {
            if (useFilt[filt])
                fprintf (w_ptr, "%6.3f ", theStar.photometry[filt] < 99. ? theStar.photometry[filt] : 99.999);
        }
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", theStar.status[0], (theStar.status[0] == 3 ? theStar.massNow[0] : 0.0), 0, theStar.wdLogTeff[0], ltau[0]);

        tempU = theStar.U;              // save so we can output the total photometry below
        if (theStar.status[0] != NSBH && theStar.status[0] != WD)
        {
            if (genrand_res53 () < fractionBinary)
            {                           // create binaries among appropriate fraction
                do
                {
                    theStar.massRatio = gen_norm (1, 0.5);
                } while (theStar.massRatio < 0 || theStar.massRatio > 1);
                theStar.U = tempU * theStar.massRatio;
                nStars++;
                massTotal += theStar.U;
                theStar.massRatio = 0.0;
            }
            else
                theStar.U = 0.0;
        }
        else
            theStar.U = 0.0;

        evolve (&theCluster, &theStar, 0);      // Evolve secondary star by itself

        fprintf (w_ptr, "%7.3f ", getMass1 (&theStar, &theCluster));    // output secondary star data
        for (filt = 0; filt < FILTS; filt++)
        {
            if (useFilt[filt])
                fprintf (w_ptr, "%6.3f ", theStar.photometry[filt] < 99. ? theStar.photometry[filt] : 99.999);
        }
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", theStar.status[0], (theStar.status[0] == 3 ? theStar.massNow[0] : 0.0), 0, theStar.wdLogTeff[0], ltau[0]);

        theStar.massRatio = theStar.U / tempU;
        theStar.U = tempU;

        evolve (&theCluster, &theStar, 0);      // Find the photometry for the whole system
        for (cmpnt = 0; cmpnt < 2; cmpnt++)
            updateCount (&theStar, cmpnt);

        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                fprintf (w_ptr, "%6.3f ", theStar.photometry[filt] < 99. ? theStar.photometry[filt] : 99.999);  // output photometry for whole system
        fprintf (w_ptr, "\n");

        // Update min and max
        if (theStar.photometry[2] < 97 && theStar.status[0] != WD && theStar.status[1] != WD)
        {
            if (theStar.photometry[2] < minV)
                minV = theStar.photometry[2];
            if (theStar.photometry[2] > maxV)
                maxV = theStar.photometry[2];
        }
    }

    ///////////////////////////////
    ///// Create brown dwarfs /////
    ///////////////////////////////

    for (i = 0; i < nBrownDwarfs; i++)
    {                           // for all systems in the cluster
        theStar.U = 0.0995 * genrand_res53 () + 0.0005; // create single stars in arbitrary mass order
        nStars++;                       // keep track of the number of stars created
        massTotal += theStar.U;
        theStar.massRatio = 0.0;
        theStar.status[0] = BD;

        evolve (&theCluster, &theStar, 0);      // given inputs, derive mags for first component

        fprintf (w_ptr, "%4d %7.4f ", i + 10001, getMass1 (&theStar, &theCluster));     // output primary star data
        for (filt = 0; filt < FILTS; filt++)
        {
            if (useFilt[filt])
                fprintf (w_ptr, "%6.3f ", theStar.photometry[filt] < 99. ? theStar.photometry[filt] : 99.999);
        }
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", theStar.status[0], 0.0, 0, theStar.wdLogTeff[0], ltau[0]);

        fprintf (w_ptr, "%7.3f ", 0.00);        // output secondary star data
        for (filt = 0; filt < FILTS; filt++)
        {
            if (useFilt[filt])
                fprintf (w_ptr, "%6.3f ", 99.999);
        }
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", DNE, 0.0, 0, 0.0, 0.0);

        for (cmpnt = 0; cmpnt < 2; cmpnt++)
            updateCount (&theStar, cmpnt);

        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                fprintf (w_ptr, "%6.3f ", theStar.photometry[filt] < 99. ? theStar.photometry[filt] : 99.999);  // output photometry for whole system
        fprintf (w_ptr, "\n");
    }

    printf ("\n Properties for cluster:");
    printf ("\n logClusAge     = %6.3f\n [Fe/H]         = %5.2f\n Y              = %5.2f", theCluster.parameter[AGE], theCluster.parameter[FEH], theCluster.parameter[YYY]);
    printf ("\n modulus        = %5.2f\n Av             = %5.2f\n WDMassUp       = %4.1f", theCluster.parameter[MOD], theCluster.parameter[ABS], theCluster.M_wd_up);
    printf ("\n fractionBinary = %5.2f\n", fractionBinary);

    printf ("\n Totals:");
    printf ("\n nSystems       = %d\n nStars         = %d\n nMSRG          = %d\n nWD            = %d", theCluster.nStars, nStars, nMSRG, nWD);
    printf ("\n nNSBH          = %d\n massTotal      = %.2f\n MSRGMassTotal  = %.2f\n wdMassTotal    = %.2f\n", nNSBH, massTotal, MSRGMassTotal, wdMassTotal);
    printf ("\n nFieldStars    = %d\n", nFieldStars * 2);


    //////////////////////////////
    ///// Create field stars /////
    //////////////////////////////

    double minFeH, maxFeH, minAge, maxAge;

    theCluster.nStars = nFieldStars;
    tempMod = theCluster.parameter[MOD];
//   tempAbs = theCluster.parameter[ABS];
    if (theCluster.evoModels.mainSequenceEvol == YALE)
        theCluster.evoModels.mainSequenceEvol = DSED;

//     setModels(&theCluster, DARTMOUTH); // Because the Yale models have trouble making field stars
    if (theCluster.evoModels.mainSequenceEvol == DSED)
    {
        minFeH = -2.5;
        maxFeH = 0.56;
        minAge = 8.4;
        maxAge = 10.17;
    }
    else
    {
        minFeH = -1.5;
        maxFeH = 0.2;
        minAge = 8.0;
        maxAge = 9.7;
    }

    i = 0;
    while (i < theCluster.nStars)
    {
        do
        {
            do
            {
                theStar.U = drawFromIMF ();
            } while (theStar.U < minMass);
            theStar.massRatio = genrand_res53 ();

            // Draw a new age and metallicity
            theCluster.parameter[AGE] = minAge + (maxAge - minAge) * genrand_res53 ();
            theCluster.parameter[FEH] = minFeH + (maxFeH - minFeH) * genrand_res53 ();

            // Determine a new distance, weighted so
            // there are more stars behind than in front
            theCluster.parameter[MOD] = tempMod - 12.0 + log10 (pow (10, (pow (pow (26.0, 3.0) * genrand_res53 (), 1.0 / 3.0))));
            //printf("%f\n",log10(pow(10,(pow(pow(5000.0,1.0)*genrand_res53(),1.0/1.0)))));


            // Calculate absorption based on rough galactic approximation and the distance.
            //theCluster.parameter[ABS] = 4*genrand_res53()*pow(10,theCluster.parameter[MOD]/5.0-2);
            //if(theCluster.parameter[ABS] < 0.0) theCluster.parameter[ABS] = 0.0;

            evolve (&theCluster, &theStar, 0);

        } while (theStar.photometry[2] < minV || theStar.photometry[2] > maxV || theStar.photometry[1] - theStar.photometry[2] < -0.5 || theStar.photometry[1] - theStar.photometry[2] > 1.7);

        fprintf (w_ptr, "%4d %7.3f ", i + 20001, getMass1 (&theStar, &theCluster));
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                fprintf (w_ptr, "%6.3f ", 99.999);
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", theStar.status[0], (theStar.status[0] == 3 ? theStar.massNow[0] : 0.0), 0, theStar.wdLogTeff[0], ltau[0]);
        fprintf (w_ptr, "%7.3f ", getMass2 (&theStar, &theCluster));
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                fprintf (w_ptr, "%6.3f ", 99.999);
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", theStar.status[1], 0.0, 0, theStar.wdLogTeff[1], ltau[1]);
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                fprintf (w_ptr, "%9.6f ", theStar.photometry[filt]);
        fprintf (w_ptr, "\n");

        i++;
    }

    fclose (w_ptr);

    return (0);
}

void updateCount (struct star *pStar, int cmpnt)
{

    switch (pStar->status[cmpnt])
    {
        case MSRG:
            nMSRG++;
            MSRGMassTotal += pStar->massNow[cmpnt];
            break;
        case WD:
            nWD++;
            wdMassTotal += pStar->massNow[cmpnt];
            break;
        case NSBH:
            nNSBH++;
            break;
    }
}
