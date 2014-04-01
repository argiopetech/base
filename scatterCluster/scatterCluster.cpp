#include <iostream>
#include <random>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <boost/format.hpp>

#include "Model.hpp"
#include "structures.hpp"
#include "Settings.hpp"
#include "MsFilterSet.hpp"

const int SCATTER_FILTS = 8;

using std::vector;
using std::cerr;
using std::cout;
using std::endl;

static double mass1, mass2, phot[SCATTER_FILTS], exptime[SCATTER_FILTS], sigma[SCATTER_FILTS];
static int stage1, stage2, starID;

double signalToNoise (double mag, double exptime, int filt);
int readLine (FILE * filePtr);
int magCutoff (int firstFilt, double brightLimit, double faintLimit);
int stageCutoff (int isFS);
int scatterPhot (double limitSigToNoise, std::mt19937);
double gen_norm (double mean, double std_dev);  // Returns a normal rv
int outputScatter (FILE * w_ptr, int isFS, double clusterMemberPrior);

vector<int> filters;

int *filterPriorMin;
int *filterPriorMax;

int main (int argc, char *argv[])
{

    int count, nr, nStars, filterSet, wdCount, i, firstFilt, nFieldStars, isFS, isBD;
    double limitSigToNoise, brightLimit, faintLimit, clusterMemberPrior;
    char filename[100], line[1000], aFilterName[10];
    FILE *r_ptr, *w_ptr;

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

    Model evoModels = makeModel(settings);

    strcpy (filename, settings.files.output.c_str());
    strcat (filename, ".sim.out");
    if ((r_ptr = fopen (filename, "r")) == NULL)
    {
        cerr << "\nFile " << filename << " was not found - exiting" << endl;
        exit (1);
    }

    //Scan header line to figure out which photometry set is being used
    //for(i=0;i<29;i++) fscanf(r_ptr,"%*s ");
    //fscanf(r_ptr,"%s",aFilterName);

    for (i = 0; i < 2; i++)
        fscanf (r_ptr, "%*s ");
    fscanf (r_ptr, "%s", aFilterName);
    i = 0;
    while (aFilterName[i] != '1')
        i++;
    aFilterName[i] = '\0';


    for (filterSet = 0; filterSet < 3; filterSet++)
    {
        if (evoModels.filterSet->getFilterName(0) == aFilterName)
            break;
    }

    fgets (line, 1000, r_ptr);  // remove rest of header line

    std::copy(settings.scatterCluster.exposures.begin(), settings.scatterCluster.exposures.end(), exptime);

    nStars = settings.simCluster.nStars;
    brightLimit = settings.scatterCluster.brightLimit;
    faintLimit = settings.scatterCluster.faintLimit;
    firstFilt = settings.scatterCluster.relevantFilt;

    limitSigToNoise = settings.scatterCluster.limitS2N;

    nFieldStars = settings.simCluster.nFieldStars;
    if (nFieldStars < 0)
        nFieldStars = 0;

    strcpy (filename, settings.files.output.c_str());
    strcat (filename, ".sim.scatter");
    if ((w_ptr = fopen (filename, "w")) == NULL)
    {
        cerr << "\nFile " << filename << " not available for writing - exiting" << endl;
        exit (1);
    }

    clusterMemberPrior = (double) nStars / (nStars + nFieldStars);
    if(clusterMemberPrior > 0.99) clusterMemberPrior = 0.99;


    std::mt19937 gen(uint32_t(settings.seed * uint32_t(2654435761))); // Applies Knuth's multiplicative hash for obfuscation (TAOCP Vol. 3)

    //Output headers
    fprintf (w_ptr, "id ");

    for (int filt = 0; filt < SCATTER_FILTS; filt++)
        if (exptime[filt] > EPS)
            fprintf (w_ptr, "%6s ", evoModels.filterSet->getFilterName (filt).c_str());

    for (int filt = 0; filt < SCATTER_FILTS; filt++)
        if (exptime[filt] > EPS)
            fprintf (w_ptr, "sig%-5s ", evoModels.filterSet->getFilterName (filt).c_str());

    fprintf (w_ptr, "mass1 massRatio stage1 CMprior useDuringBurnIn\n");

    count = 0;
    wdCount = 0;

    isFS = 0;
    isBD = 0;

    //Fix this line to read in however many bands there are
    while ((nr = readLine (r_ptr)) != EOF)
    {
        // simCluster creates brown dwarfs with #'s starting at 10,001
        if (starID > 10000 && !isBD)
        {
            if (count < nStars)
            {
                cerr << " Warning - not enough (" << count << ") non-RGB stars in input file to keep desired number (" << nStars << ") of stars." << endl;
            }
            count = 0;
            nStars = 100000;
            isBD = 1;
        }

        // simCluster creates field stars with #'s starting at 20,001
        if (starID > 20000 && !isFS)
        {
            count = 0;
            isFS = 1;
        }

        if (magCutoff (firstFilt, brightLimit, faintLimit) == 0)
            continue;
        if (stageCutoff (isFS) == 0)
            continue;
        if (scatterPhot (limitSigToNoise, gen) == 0)
            continue;
        if (mass1 < 0.25 && mass2 < 0.25 && stage1 != BD)
            continue;                   // limit on modelSet=4 is 0.4 Mo - what to do?
        //if(stage1 == BD) continue;
        wdCount += outputScatter (w_ptr, isFS, clusterMemberPrior);
        count++;

        // If we have enough stars of this type
        if (count == nStars)
        {
            // If these are field stars, we're done
            if (isFS)
                break;
            // If not...
            else
            {
                // Skip lines in the file until you get to...
                while ((nr = fscanf (r_ptr, "%d ", &starID)) != EOF)
                {
                    //...the brown dwarfs (which start at 10001)...
                    if (starID == 10001)
                    {
                        isBD = 1;
                        break;
                    }
                    //...or the field stars (which start at 20001)
                    if (i == 20001)
                    {
                        isFS = 1;
                        isBD = 1;
                        break;
                    }
                    fgets (line, 1000, r_ptr);
                }
                fseek (r_ptr, -7, SEEK_CUR);

                // If this is the first field star, reset the nStars to
                // nFieldStars and read in nFieldStars in the main loop
                if (isFS)
                {
                    nStars = nFieldStars;
                    count = 0;
                    continue;
                }
                else if (isBD)
                {
                    nStars = 100000000;
                    count = 0;
                    continue;
                }
            }
        }
    }
    if (count < nFieldStars)
    {
        cerr << " Warning - not enough (" << count << ") field stars in input file to keep desired number (" << nFieldStars << ") of stars." << endl;
    }

    cerr << " There " << (wdCount == 1 ? "is " : "are " )<< wdCount << " WD%s in this scatter file, " << filename << endl;

    fclose (r_ptr);
    fclose (w_ptr);

    return (0);
}


static double s2nCoeffs[][2] = {
    {9.33989, 0.3375778},       // U
    {10.0478, 0.3462758},       // B
    {10.48098, 0.368201},       // V
    {10.71151, 0.3837847},      // R
    {10.61035, 0.3930941},      // I
    {9.282385, 0.386258},       // J
    {9.197463, 0.3970419},      // H
    {9.024068, 0.3985604},      // K
    {9.024068, 0.3985604},      // IRAC Blue
    {9.024068, 0.3985604},      // IRAC Red
    {9.024068, 0.3985604},      // Band 1
    {9.024068, 0.3985604},      // Band 2
    {9.024068, 0.3985604},      // Band 3
    {9.024068, 0.3985604}       // Band 4
};

double signalToNoise (double mag, double exptime, int filter)
/*
  This is an approximation to the results one would obtain in one hour with the KPNO
  4m + Mosaic (UBVRI) or Flamingos (JHK) per band, assuming dark time, seeing=1.1",
  airmass=1.2, and then scaling from their by sqrt(exptime).  I further approximated
  the CCDTIME results (run on the NOAO webste) with linear fits of mag vs. log(S/N).
*/
{

    double s2n, logS2N;

    if (filter >= SCATTER_FILTS || filter < 0)
    {
        cerr << "Filter (" << filter << ") out of range - exiting" << endl;
        exit (1);
    }

    // Scatter BD photometry at 5%
    if (filter > 7)
        return 1.0 / 0.05;

    logS2N = s2nCoeffs[filter][0] - s2nCoeffs[filter][1] * mag;
    s2n = pow (10., logS2N);
    s2n *= sqrt (exptime);

    return (s2n);

}

int readLine (FILE * filePtr)
{

    int nr, filt;

    fscanf (filePtr, "%d %lf ", &starID, &mass1);

    for (filt = 0; filt < SCATTER_FILTS; filt++)
    {
        fscanf (filePtr, "%*lf ");
    }

    fscanf (filePtr, "%d %*f %*d %*f %*f %lf ", &stage1, &mass2);
    for (filt = 0; filt < SCATTER_FILTS; filt++)
        fscanf (filePtr, "%*f ");
    nr = fscanf (filePtr, "%d %*f %*d %*f %*f ", &stage2);
    for (filt = 0; filt < SCATTER_FILTS; filt++)
    {
        nr = fscanf (filePtr, "%lf ", &phot[filt]);
    }

    if (starID == EOF)
    {
        cerr << "\nFatal error in readLine.  Exiting." << endl;
        exit (1);
    }

    if (nr == EOF)
        return nr;
    return starID;

}

int magCutoff (int firstFilt, double brightLimit, double faintLimit)
{
    if (phot[firstFilt] > 99. && stage1 != BD)
        return 0;                       // check if real object.  if not, skip
    if (phot[firstFilt] < brightLimit && stage1 != BD)
        return 0;                       // typically used to remove red giants
    if (phot[firstFilt] > faintLimit && stage1 != BD)
        return 0;                       // typically used to remove WDs and lower MS
    return 1;
}

int stageCutoff (int isFS)
{
    if (stage1 == NSBH || stage2 == NSBH)
        return 0;
    if (stage1 == WD && mass2 > 0.0 && !isFS)
        return 0;                       // TEMPORARY KLUDGE -- ignore binaries of MS/RG + WDs and WD + WD
    if (stage2 == WD && mass1 > 0.0 && !isFS)
        return 0;
    return 1;
}

int scatterPhot (double limitSigToNoise, std::mt19937 gen)
{
    int filt;
    double sigToNoise;

    for (filt = (stage1 == BD ? 8 : 0); filt < (stage1 == BD ? SCATTER_FILTS : 8); filt++)
    {
        if (exptime[filt] < EPS)
            continue;
        sigToNoise = signalToNoise (phot[filt], exptime[filt], filt);
        if (sigToNoise < limitSigToNoise)
        {                               // large photometric errors can lock mcmc during burnin
            cerr << boost::format("Warning: star %4d, mass1=%.3f, stage1=%d, filter=%d, S/N = %.3f < %.1f (user limit) - skipping.") % starID % mass1 % stage1 % filt % sigToNoise % limitSigToNoise << endl;
            return 0;
        }

        sigma[filt] = 1. / (sigToNoise);
        if (sigma[filt] < 0.01)
            sigma[filt] = 0.01;

        std::normal_distribution<double> normDist(0., sigma[filt]);
        phot[filt] += normDist(gen);
    }
    return 1;
}

int outputScatter (FILE * w_ptr, int isFS, double clusterMemberPrior)
{
    int tempStage, filt;
    double tempMass, tempMassRatio;

    if (mass1 > mass2)
    {                           // use to set starter mass for mcmc
        tempMass = mass1;               // (since higher mass star dominates the photometry)
        tempMassRatio = mass2 / mass1;
        tempStage = stage1;
    }
    else
    {
        tempMass = mass2;
        tempMassRatio = mass1 / mass2;
        tempStage = stage2;
    }

    //Fix outputs
    fprintf (w_ptr, "%6d ", starID);
    for (filt = 0; filt < SCATTER_FILTS; filt++)
        if (exptime[filt] > EPS)
            fprintf (w_ptr, "%6.3f ", phot[filt]);
    for (filt = 0; filt < 8; filt++)
    {
        if (exptime[filt] > EPS)
        {
            if (stage1 == BD)
                fprintf (w_ptr, "%8.5f ", -1.0);
            else
                fprintf(w_ptr,"%8.6f ",sigma[filt]);
        }
    }
    for (filt = 8; filt < SCATTER_FILTS; filt++)
    {
        if (exptime[filt] > EPS)
        {
            if(stage1 == BD)
                fprintf(w_ptr,"%8.6f ",sigma[filt]);
            else
                fprintf (w_ptr, "%8.5f ", -1.0);
        }
    }

    fprintf (w_ptr, "%8.3f %6.3f %3d %6.3f   %d\n", tempMass, (tempMassRatio > 0.001 ? tempMassRatio : 0.00), tempStage, clusterMemberPrior, !isFS);

    if (tempStage == WD && !isFS)
        return 1;
    else
        return 0;
}
