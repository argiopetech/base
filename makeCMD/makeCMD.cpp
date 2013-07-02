#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <cmath>
#include <cstring>
#include <unistd.h>

#include "evolve.hpp"
#include "structures.hpp"
#include "loadModels.hpp"
#include "Settings.hpp"
#include "ifmr.hpp"
#include "FilterSet.hpp"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

const double  FS_NUM       = 0.0;    // mcmc outputs negative mass for field star
const int  COL_MAX         = 1000;    // max number of cluster stars
const int  CLUS_READ       = 9;
const int  CLUS_STAT_WRITE = 10;
const int  MASS1_FILE      = 11;
const int  MASS2_FILE      = 12;
const int  CMD_FILE        = 13;
const int  DEBUG_FILE      = 16;

// Used by evolve.c
double ltau[2], wdLogTeff[2];
int aFilt = 0;

// Used by a bunch of different functions.
int needMassNow = 1, useFilt[FILTS];   //, numFilts;

void openOutputFiles (FILE ** filePtr, char *filename, int fileType);

int main (int argc, char *argv[])
{

    /////////////////////////////////
    /////// Declare Variables ///////
    /////////////////////////////////

    int j, nr, p, m, row, filt, isClusterMember[2][COL_MAX], cm[COL_MAX], stuckMass[2][COL_MAX], stuckParam[NPARAMS], tempCols, useParam[NPARAMS], numStepSizeMass[2][COL_MAX], iMag;
    double mass[2][COL_MAX], prevMass[2][COL_MAX], mean = 0.0, minMag, maxMag, sumMass[2][COL_MAX], sumSquaresMass[2][COL_MAX], minMass[2][COL_MAX], maxMass[2][COL_MAX], thisStep, meanStepSizeMass[2][COL_MAX], meanMass[2][COL_MAX], variance = 0.0, sigma, N, csProb;
    double param[NPARAMS], meanParam[NPARAMS], varParam[NPARAMS], minParam[NPARAMS], maxParam[NPARAMS], prevParam[NPARAMS], sumStepSizeParam[NPARAMS];
    double sumPhot[COL_MAX][FILTS], sumSquaresPhot[COL_MAX][FILTS];
    char line[240], filename[100], **starName;

    const char *paramNames[] = { "logAge    ",
                                 "Y         ",
                                 "[Fe/H]    ",
                                 "modulus   ",
                                 "absorption"
    };

    FILE *rMassPtr[2], *rClusterPtr, *rDataPtr;
    FILE *wClusterStatPtr, *wCmdPtr, *wDebugPtr;

    Cluster theCluster;
    vector<Star> stars;

    Settings settings;

    settings.fromCLI (argc, argv);
    if (!settings.files.config.empty())
    {
        settings.fromYaml (settings.files.config);
    }
    else
    {
        settings.fromYaml ("base9.yaml");
    }

    settings.fromCLI (argc, argv);

    const Model evoModels = makeModel(settings);

    ////////////////////////////////////////////
    /////// Open files to read and write ///////
    ////////////////////////////////////////////

    strcpy (filename, settings.files.scatter.c_str());
    if ((rDataPtr = fopen (filename, "r")) == NULL)
    {
        cerr << "***Error: file " << filename << " was not found.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    minMag = settings.cluster.minMag;
    maxMag = settings.cluster.maxMag;
    iMag = settings.cluster.index;

    if (iMag < 0 || iMag > 2)
    {
        cerr << "***Error: " << iMag << " not a valid magnitude index.  Choose 0,1,or 2.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    // This is a leftover from Base8 and may not work with the current cluster files
    strcpy (filename, settings.files.output.c_str());
    openOutputFiles (&rClusterPtr, filename, CLUS_READ);
    openOutputFiles (&wClusterStatPtr, filename, CLUS_STAT_WRITE);
    openOutputFiles (&wCmdPtr, filename, CMD_FILE);
    openOutputFiles (&wDebugPtr, filename, DEBUG_FILE);


    ///////////////////////////////////////
    //////// Determine how many stars /////
    ////////// and allocate memory ////////
    ///////////////////////////////////////

    for (m = 0; m < 2; m++)
    {
        openOutputFiles (&rMassPtr[m], filename, MASS1_FILE + m);
        fscanf (rMassPtr[m], "%*s %*s %*s %*s %*s %*s %d", ((m) ? &(theCluster.nStars) : &tempCols));   // Read in the number of stars
    }

    if (theCluster.nStars < 1)
    {
        cerr << "\nNeed at least a one column file - exiting." << endl;
        exit (1);
    }
    if (theCluster.nStars > COL_MAX)
    {
        cerr << "\nExceeded the limit of the number of cluster stars in a file - exiting." << endl;
        exit (1);
    }

    if (theCluster.nStars != tempCols)
    {
        cerr << "\nDifferent numbers of stars in .mass1 and .mass2 files" << endl;
        exit (1);
    }

    stars.resize(theCluster.nStars);

    starName = new char*[theCluster.nStars]();

    for (j = 0; j < theCluster.nStars; j++)
    {
        starName[j] = new char[100]();
    }

    ////////////////////////////////////////
    //// Determine which parameters are ////
    //// included in the .cluster file /////
    ////////////////////////////////////////

    // Skip first line of .cluster file
    fgets (line, 240, rClusterPtr);

    //Read in second line
    fscanf (rClusterPtr, "%*s ");
    for (p = 0; p < NPARAMS; p++)
    {
        fscanf (rClusterPtr, "%d %lf ", &(useParam[p]), &(param[p]));
    }
    fgets (line, 240, rClusterPtr);


    ////////////////////////////////////////////
    ////// Determine which filters to use //////
    ////// read header from .scatter file //////
    ////////////////////////////////////////////


    fgets (line, 300, rDataPtr);        // skip first header line

    for (filt = 0; filt < FILTS; filt++)
    {
        fscanf (rDataPtr, "%d ", &useFilt[filt]);
        if (useFilt[filt])
        {
            const_cast<Model&>(evoModels).numFilts++;
            aFilt = filt;               // Sets this to a band we know we are using (for evolve)
        }
    }
    fgets (line, 240, rDataPtr);        // and next header line


    ///////////////////////////////////
    ///// Read in other variables /////
    ///////// and load models /////////
    ///////////////////////////////////

    theCluster.M_wd_up = settings.whiteDwarf.M_wd_up;
    // verbose = settings.verbose;

    // if (verbose < 0 || verbose > 2)
    //     verbose = 1;            // give standard feedback if incorrectly specified

    loadModels (&theCluster, evoModels, settings);      // read in stellar evol & WD models


    //////////////////////////////
    /////// Output headers ///////
    //////////////////////////////

    fprintf (wCmdPtr, "     Star    ");
    for (m = 1; m <= 2; m++)
    {
        fprintf (wCmdPtr, "meanmass%d  sigmass%d   minmass%d   maxmass%d stuck%d meanStep%d ", m, m, m, m, m, m);
        if (m == 1)
            fprintf (wCmdPtr, "meanWDMass%d sigWDmass%d ", m, m);
        fprintf (wCmdPtr, "CSprob%d Niter%d   ", m, m);
    }
    for (filt = 0; filt < FILTS; filt++)
    {
        if (useFilt[filt])
        {
            char tmpString[100] = "mean\0";

            strcat (tmpString, getFilterName (filt));
            fprintf (wCmdPtr, "%8s ", tmpString);
            tmpString[0] = '\0';
            strcat (tmpString, "var");
            strcat (tmpString, getFilterName (filt));
            fprintf (wCmdPtr, "  %8s   ", tmpString);
        }
    }
    fprintf (wCmdPtr, "\n");


    ////////////////////////////////
    ///// Initialize variables /////
    ////////////////////////////////

//   m=0;
    nr = 0;
    row = 0;
    for (p = 0; p < NPARAMS; p++)
    {
        minParam[p] = HUGE_VAL;
        maxParam[p] = -HUGE_VAL;
        stuckParam[p] = 0;
        meanParam[p] = 0;
        varParam[p] = 0;
        sumStepSizeParam[p] = 0;
    }
    for (j = 0; j < theCluster.nStars; j++)
    {
        for (filt = 0; filt < FILTS; filt++)
        {
            sumPhot[j][filt] = 0.0;
            sumSquaresPhot[j][filt] = 0.0;
        }
    }

    for (m = 0; m < 2; m++)
    {
        for (j = 0; j < COL_MAX; j++)
        {                               // initialize variables
            sumMass[m][j] = 0.0;
            sumSquaresMass[m][j] = 0.0;
            isClusterMember[m][j] = 0;
            minMass[m][j] = 100;        // impossibly high min mcmc value for each mass
            maxMass[m][j] = -100;       // impossibly low max mcmc value for each mass
            prevMass[m][j] = -100.;     // impossible value, should never be used
            stuckMass[m][j] = 0;
            meanStepSizeMass[m][j] = 0.;
            numStepSizeMass[m][j] = 0;
        }
    }


    ////////////////////////////
    ///////// Main loop ////////
    ////////////////////////////

    while (1)
    {
        // Read a line in from both mass files and the cluster file
        for (j = 0; j < theCluster.nStars; j++)
        {
            if ((nr = fscanf (rMassPtr[0], "%lf", &mass[0][j])) == EOF)
                break;
            fscanf (rMassPtr[1], "%lf", &mass[1][j]);
            //And the .scatter file
            do
            {
                fscanf (rDataPtr, "%s ", starName[j]);
                for (p = 0; p < evoModels.numFilts; p++)
                    fscanf (rDataPtr, "%lf ", &(stars.at(j).obsPhot[p]));
                for (p = 0; p < evoModels.numFilts; p++)
                    fscanf (rDataPtr, "%*f ");
                fscanf (rDataPtr, "%*f %*f %d ", &stars.at(j).status[0]);
                fgets (line, 240, rDataPtr);
            } while (stars.at(j).status[0] == MSRG && (stars.at(j).obsPhot[iMag] < minMag || stars.at(j).obsPhot[iMag] > maxMag));
        }
        if (nr == EOF)
            break;

        fscanf (rClusterPtr, "%*d ");
        // Store stats for each cluster parameter
        for (p = 0; p < NPARAMS; p++)
        {
            if (useParam[p])
                fscanf (rClusterPtr, "%lf ", &param[p]);
            meanParam[p] = (row * meanParam[p] + param[p]) / (row + 1); // Calculate mean and variance recursively
            if (row > 0)
                varParam[p] = row * varParam[p] / (row + 1) + pow ((param[p] - meanParam[p]), 2) / row;
            if (param[p] > maxParam[p])
                maxParam[p] = param[p];
            if (param[p] < minParam[p])
                minParam[p] = param[p];

            thisStep = fabs (param[p] - prevParam[p]);  // accumulate correlation and sticking stats
            if (thisStep < EPS)
                stuckParam[p]++;
            sumStepSizeParam[p] += thisStep;
            prevParam[p] = param[p];
            theCluster.parameter[p] = param[p];
        }


        // Store stats for each mass for each star
        for (j = 0; j < theCluster.nStars; j++)
        {
            if (mass[0][j] > EPS)
                cm[j] = 1;
            else
                cm[j] = 0;
            if (cm[j])
            {                           // if it is a cluster star
                for (m = 0; m < 2; m++)
                {                       // For each component
                    sumMass[m][j] += mass[m][j];        // for mean, sigma, N
                    sumSquaresMass[m][j] += mass[m][j] * mass[m][j];

                    if (minMass[m][j] > mass[m][j])
                        minMass[m][j] = mass[m][j];     // min & max
                    if (maxMass[m][j] < mass[m][j])
                        maxMass[m][j] = mass[m][j];

                    if (prevMass[m][j] > 0)
                    {
                        thisStep = fabs (mass[m][j] - prevMass[m][j]);  // accumulate correlation and sticking stats
                        if (thisStep < EPS)
                            stuckMass[m][j]++;
                        meanStepSizeMass[m][j] += thisStep;
                        numStepSizeMass[m][j]++;
                    }

                    prevMass[m][j] = mass[m][j];
                    isClusterMember[m][j]++;    // Keep track of how many iterations we're averaging for each star

                    stars.at(j).U = mass[0][j];
                    stars.at(j).massRatio = mass[1][j] / mass[0][j];
                }
            }
        }

        for (auto s : stars)
            evolve (theCluster, evoModels, s);

        for (j = 0; j < theCluster.nStars; j++)
        {
            if (cm[j])
            {                           // if it is a cluster star
                for (filt = 0; filt < evoModels.numFilts; filt++)
                {
                    sumPhot[j][filt] += stars.at(j).photometry[filt];
                    sumSquaresPhot[j][filt] += stars.at(j).photometry[filt] * stars.at(j).photometry[filt];
                }
            }
        }
        row++;
    }

    /////////////////////////////////////////////
    //////// Calculate stats and output /////////
    /////////////////////////////////////////////

    // Calculate mass stats for each star, for each component
    for (j = 0; j < theCluster.nStars; j++)
    {                           // calc mean, sigma values for each star

        fprintf (wCmdPtr, "%9s  ", starName[j]);
        for (m = 0; m < 2; m++)
        {
            N = (double) isClusterMember[m][j];
            if (N > 1.0)
            {
                variance = (N * sumSquaresMass[m][j] - sumMass[m][j] * sumMass[m][j]) / (N * (N - 1));
                sigma = sqrt (fabs (variance));
            }
            else
                sigma = 0.0;
            if (N >= 0.5)
            {
                meanMass[m][j] = sumMass[m][j] / N;     // else mean = 0, anyway
                if (numStepSizeMass[m][j] > 0)
                    meanStepSizeMass[m][j] /= numStepSizeMass[m][j];
            }

            csProb = N / row;           // row = mcmc iterations, N = # of times it was a cluster member

            // Output mass stats
            fprintf (wCmdPtr, "%11.5f %9.5f %10.5f %10.5f %6d  %8.6f  ", meanMass[m][j], sigma, minMass[m][j], maxMass[m][j], stuckMass[m][j], meanStepSizeMass[m][j]);
            if (m == 0)
            {
                if (stars.at(j).status[0] == WD)
                    // fprintf(wCmdPtr,"%10.5f  %9.3e  ", intlFinalMassReln(meanMass[m][j], theCluster.evoModels.IFMR),
                    //         intlFinalMassReln(meanMass[m][j], theCluster.evoModels.IFMR) -
                    //         intlFinalMassReln(meanMass[m][j] - sigma, theCluster.evoModels.IFMR));
                    fprintf (wCmdPtr, "%10.5f  %9.3e  ", intlFinalMassReln (theCluster, evoModels, meanMass[m][j]), intlFinalMassReln (theCluster, evoModels, meanMass[m][j]) - intlFinalMassReln (theCluster, evoModels, meanMass[m][j] - sigma));
                else
                    fprintf (wCmdPtr, "%10.5f  %9.3e  ", 0.0, 0.0);
            }
            fprintf (wCmdPtr, "%6.4f %6.0f ", csProb, N);


            // Calculate photometry stats
            if (m == 1)
            {
                for (filt = 0; filt < evoModels.numFilts; filt++)
                {
                    if (N > 1.0)
                    {
                        variance = (N * sumSquaresPhot[j][filt] - sumPhot[j][filt] * sumPhot[j][filt]) / (N * (N - 1));
                    }
                    if (N >= 0.5)
                    {
                        mean = sumPhot[j][filt] / N;
                    }
                    else
                        mean = 99.999;

                    fprintf (wCmdPtr, "%10f %10.4e ", mean < 99 ? mean : 99.999, mean < 99 ? variance : 0.0);
                }
                fprintf (wCmdPtr, "\n");
            }
        }
    }

    // Output stats for cluster parameters
    fprintf (wClusterStatPtr, "Parameter     mean     sigma       min      max    stuck  meanStep\n");
    for (p = 0; p < NPARAMS; p++)
        fprintf (wClusterStatPtr, "%s %8.5f %9.5f %10.5f %8.5f %6d  %8.5f\n", paramNames[p], meanParam[p], sqrt (varParam[p]), minParam[p], maxParam[p], stuckParam[p], sumStepSizeParam[p] / (row - 1));

    // Create a simulated cluster based on the cluster stats and output for comparison
    theCluster.nStars = (int) (theCluster.M_wd_up * 100) + 1;
    for (p = 0; p < NPARAMS; p++)
        theCluster.parameter[p] = meanParam[p];

    stars.resize(theCluster.nStars);

    for (j = 0; j < theCluster.nStars; j++)
    {
        stars.at(j).U = j * .01;
        stars.at(j).massRatio = 0.0;
    }

    for (auto s : stars)
        evolve (theCluster, evoModels, s);

    fprintf (wDebugPtr, " mass stage1");
    for (filt = 0; filt < FILTS; filt++)
        if (useFilt[filt])
            fprintf (wDebugPtr, "          %s", getFilterName (filt));
    fprintf (wDebugPtr, "\n");

    double prevV = 100;

    for (j = 0; j < theCluster.nStars; j++)
    {
        if (stars.at(j).photometry[0] < 90)
        {
            if (stars.at(j).status[0] == MSRG)
            {
                if (prevV > stars.at(j).photometry[0])
                {
                    fprintf (wDebugPtr, "%5.2f %6d ", stars.at(j).U, stars.at(j).status[0]);
                    for (filt = 0; filt < evoModels.numFilts; filt++)
                        fprintf (wDebugPtr, "%10f ", stars.at(j).photometry[filt]);
                    fprintf (wDebugPtr, "\n");
                    prevV = stars.at(j).photometry[0];
                }
                else
                {
                    prevV = 0.0;
                    continue;
                }
            }
            else
            {
                fprintf (wDebugPtr, "%5.2f %3d ", stars.at(j).U, stars.at(j).status[0]);
                for (filt = 0; filt < evoModels.numFilts; filt++)
                    fprintf (wDebugPtr, "%10f ", stars.at(j).photometry[filt]);
                fprintf (wDebugPtr, "\n");
            }
        }
    }

    fclose (wClusterStatPtr);
    fclose (rClusterPtr);
    fclose (rMassPtr[0]);
    fclose (rMassPtr[1]);
    fclose (wCmdPtr);

    return (0);
}


void openOutputFiles (FILE ** filePtr, char *filename, int fileType)
{
    char tmpfile[100];
    const char *mode = "r";

    strcpy (tmpfile, filename);
    switch (fileType)
    {
        // output files differentiated by file name extension
        case CLUS_READ:
            strcat (tmpfile, ".cluster");
            break;
        case MASS1_FILE:
            strcat (tmpfile, ".mass1");
            break;
        case MASS2_FILE:
            strcat (tmpfile, ".mass2");
            break;
        case CLUS_STAT_WRITE:
            strcat (tmpfile, ".cluster.stat");
            mode = "w";
            break;
        case CMD_FILE:
            strcat (tmpfile, ".cmd");
            mode = "w";
            break;
        case DEBUG_FILE:
            strcat (tmpfile, ".cmd.debug");
            mode = "w";
            break;
        default:
            cerr << "***Error: Bad file choice in openOutputFiles().***" << endl;
            cerr << "[Exiting...]" << endl;
            exit (0);

    }

    cout << "Reading file  : " << tmpfile << " (" << mode << ")" << endl;
    if ((*filePtr = fopen (tmpfile, mode)) == NULL)
    {
        cerr << "***Error: File " << tmpfile << " was not available. " << mode << " ***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }
}
