#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include <unistd.h>
#include <boost/format.hpp>

#include "constants.hpp"
#include "mpiMcmc.hpp"
#include "mpiFuncs.hpp"
#include "evolve.hpp"
#include "loadModels.hpp"
#include "msRgbEvol.hpp"
#include "gBergMag.hpp"
#include "wdCooling.hpp"
#include "densities.hpp"
#include "decide.hpp"
#include "samplers.hpp"
#include "leastSquares.hpp"
#include "mt19937ar.hpp"
#include "Settings.hpp"
#include "FilterSet.hpp"

using std::array;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::ofstream;
using std::isfinite;

/*** global variables ***/
/* Used by evolve.c */
//double ltau[2];
int aFilt = -1;

/* Used in densities.c. */
double filterPriorMin[FILTS];
double filterPriorMax[FILTS];
extern double ageLimit[2];      /* Defined in evolve.c, set in loadModels. */
extern struct globalIso isochrone;

int needMassNow = 0, useFilt[FILTS];


/*
 * Read data
 */
void readCmdData (Chain &mc, struct ifmrMcmcControl &ctrl, const Model &evoModels)
{
    char line[300];
    double tempSigma;
    int filt, i;
    char *pch, sig[] = "sig", comp[] = "   ";

    //Parse the header of the file to determine which filters are being used
    ctrl.rData.getline(line, 300);     // Read in the header line

    pch = strtok (line, " ");   // split the string on these delimiters into "tokens"

    while (pch != NULL)
    {
        pch = strtok (NULL, " ");       // Ignore the first token (which is "id") and move
        // to the next (which should be the first filter name)
        strncpy (comp, pch, 3); // copy the first three letters into the dummy string 'comp'
        if (strcmp (comp, sig) == 0)
            break;                      // and check to see if they are 'sig'.  If they are, there are no more filters

        for (filt = 0; filt < FILTS; filt++)
        {                               // Otherwise check to see what this filter's name is
            if (strcmp (pch, getFilterName (filt)) == 0)
            {
                ctrl.useFilt[filt] = 1;
                const_cast<Model&>(evoModels).numFilts++;
                if (aFilt < 0)
                    aFilt = filt;               // Sets this to a band we know we are using (for evolve)
                break;
            }
        }
    }

    for (i = 0; i < FILTS; i++)
    {
        if (ctrl.useFilt[i])
        {
            ctrl.numFilts++;
            if (aFilt < 0)
                aFilt = i;              // Sets this to a band we know we are using (for evolve)
        }
    }

    // This loop reads in photometry data
    // It also reads a best guess for the mass
    int j = 0;
    int moreStars = 1;          // true

    // why is this necessary???
    mc.stars.clear();

    while (moreStars)
    {
        ctrl.rData >> line;

        if (ctrl.rData.eof())
            break;

        mc.stars.emplace_back();

        for (i = 0; i < ctrl.numFilts; i++)
        {
            ctrl.rData >> mc.stars[j].obsPhot[i];

            if (mc.stars[j].obsPhot[i] < ctrl.filterPriorMin[i])
            {
                ctrl.filterPriorMin[i] = mc.stars[j].obsPhot[i];
            }

            if (mc.stars[j].obsPhot[i] > ctrl.filterPriorMax[i])
            {
                ctrl.filterPriorMax[i] = mc.stars[j].obsPhot[i];
            }
        }

        // copy to global variables
        for (i = 0; i < ctrl.numFilts; i++)
        {
            filterPriorMin[i] = ctrl.filterPriorMin[i];
            filterPriorMax[i] = ctrl.filterPriorMax[i];
        }
        for (i = 0; i < ctrl.numFilts; i++)
        {
            ctrl.rData >> tempSigma;
            mc.stars[j].variance[i] = tempSigma * fabs (tempSigma);
            // The fabs() keeps the sign of the variance the same as that input by the user for sigma
            // Negative sigma (variance) is used to signal "don't count this band for this star"
        }

        ctrl.rData >> mc.stars[j].U >> mc.stars[j].massRatio >> mc.stars[j].status[0] >> mc.stars[j].clustStarPriorDens >> mc.stars[j].useDuringBurnIn;

        if (mc.stars[j].status[0] == 3 || (mc.stars[j].obsPhot[ctrl.iMag] >= ctrl.minMag && mc.stars[j].obsPhot[ctrl.iMag] <= ctrl.maxMag))
        {
            j++;
        }
    }
    mc.clust.nStars = j;

    for (j = 0; j < mc.clust.nStars; j++)
    {
        mc.stars[j].massRatio = 0.0;
    }

    // copy to global values
    for (i = 0; i < FILTS; i++)
    {
        useFilt[i] = ctrl.useFilt[i];
    }

    assert (mc.stars.size() == mc.clust.nStars);
} /* readCmdData */



void propClustBigSteps (Cluster &clust, struct ifmrMcmcControl const &ctrl)
{
    /* DOF defined in densities.h */
    double scale = 5.0;
    int p;

    for (p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar[p] > EPSILON)
        {
            clust.parameter[p] += sampleT (scale * scale * clust.stepSize[p] * clust.stepSize[p], DOF);
        }
    }
}

void propClustIndep (Cluster &clust, struct ifmrMcmcControl const &ctrl)
{
    /* DOF defined in densities.h */
    int p;

    for (p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar[p] > EPSILON)
        {
            clust.parameter[p] += sampleT (clust.stepSize[p] * clust.stepSize[p], DOF);
        }
    }
}

void propClustCorrelated (Cluster &clust, struct ifmrMcmcControl const &ctrl)
{
    /* DOF defined in densities.h */
    double indepProps[NPARAMS] = { 0.0 };
    double corrProps[NPARAMS] = { 0.0 };

    int p, k;

    for (p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar[p] > EPSILON)
        {
            indepProps[p] = sampleT (1.0, DOF);
        }
    }
    for (p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar[p] > EPSILON)
        {
            for (k = 0; k <= p; k++)
            {                           /* propMatrix is lower diagonal */
                if (ctrl.priorVar[k] > EPSILON)
                {
                    corrProps[p] += ctrl.propMatrix[p][k] * indepProps[k];
                }
            }
            clust.parameter[p] += corrProps[p];
        }
    }
}


void printHeader (ofstream &file, array<double, NPARAMS> const &priors)
{
    const char *paramNames[] = { "    logAge",
                                 "         Y",
                                 "       FeH",
                                 "   modulus",
                                 "absorption",
                                 " IFMRconst",
                                 "   IFMRlin",
                                 "  IFMRquad"
    };

    for (int p = 0; p < NPARAMS; p++)
    {
        if (priors.at(p) > EPSILON || p == MOD || p == FEH || p == ABS)
        {
            file << paramNames[p] << ' ';
        }
    }
    file << "logPost" << endl;
}

int main (int argc, char *argv[])
{
    int accept = 0, reject = 0;
    int increment;

    double logPostCurr;
    double logPostProp;
    double fsLike;

    Chain mc;
    struct ifmrMcmcControl ctrl;
    Cluster propClust;

    array<double, 2> ltau;
    array<double, N_MS_MASS1 * N_MS_MASS_RATIO> msMass1Grid;
    array<double, N_MS_MASS1 * N_MS_MASS_RATIO> msMassRatioGrid;
    array<double, N_WD_MASS1> wdMass1Grid;

    Matrix<double, NPARAMS, nSave> params;

    Settings settings;

    // Setup settings
    {
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
    }

    const Model evoModels = makeModel(settings);

    increment = settings.mpiMcmc.burnIter / (2 * nSave);

    initIfmrMcmcControl (mc, ctrl, evoModels, settings);

    for (int p = 0; p < NPARAMS; p++)
    {
        mc.clust.priorVar[p] = ctrl.priorVar[p];
        mc.clust.priorMean[p] = ctrl.priorMean[p];
    }

    readCmdData (mc, ctrl, evoModels);

    initChain (mc, ctrl, evoModels, ltau);

    for (int i = 0; i < mc.clust.nStars; i++)
    {
        mc.stars[i].isFieldStar = 0;
        mc.stars[i].boundsFlag = 0;
    }

    initMassGrids (msMass1Grid, msMassRatioGrid, wdMass1Grid, mc);

    double logFieldStarLikelihood = 0.0;

    if (mc.clust.nStars > 1)
    {
        for (int filt = 0; filt < ctrl.numFilts; filt++)
        {
            logFieldStarLikelihood -= log (ctrl.filterPriorMax[filt] - ctrl.filterPriorMin[filt]);
        }
        fsLike = exp (logFieldStarLikelihood);
    }
    else
    {
        logFieldStarLikelihood = -HUGE_VAL;
        fsLike = 0;
    }

    cout << "Bayesian analysis of stellar evolution" << endl;

    /* set current log posterior to -HUGE_VAL */
    /* will cause random starting value */
    logPostCurr = -HUGE_VAL;

    // Run Burnin
    ctrl.burninFile.open(ctrl.clusterFilename + ".burnin");
    if (!ctrl.burninFile)
    {
        cerr << "***Error: File " << ctrl.clusterFilename << " was not available for writing.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    printHeader (ctrl.burninFile, ctrl.priorVar);

    for (int iteration = 0; iteration < ctrl.burnIter; iteration++)
    {
        propClust = mc.clust;

        if (iteration < ctrl.burnIter / 2)
        {
            propClustBigSteps (propClust, ctrl);
        }
        else
        {
            propClustIndep (propClust, ctrl);
        }

        if (ctrl.priorVar[ABS] > EPSILON)
        {
            propClust.parameter[ABS] = fabs (propClust.parameter[ABS]);
        }
        if (ctrl.priorVar[IFMR_SLOPE] > EPSILON)
        {
            propClust.parameter[IFMR_SLOPE] = fabs (propClust.parameter[IFMR_SLOPE]);
        }

        logPostProp = logPostStep (mc, evoModels, wdMass1Grid, propClust, fsLike, ltau);

        /* accept/reject */
        if (acceptClustMarg (logPostCurr, logPostProp, ltau))
        {
            mc.clust = propClust;
            logPostCurr = logPostProp;
            accept++;
        }
        else
        {
            reject++;
        }
        /* save draws to estimate covariance matrix for more efficient Metropolis */
        if (iteration >= ctrl.burnIter / 2 && iteration < ctrl.burnIter)
        {
            if (iteration % increment == 0)
            {
                /* save draws */
                for (int p = 0; p < NPARAMS; p++)
                {
                    if (ctrl.priorVar[p] > EPSILON)
                    {
                        params.at(p).at((iteration - ctrl.burnIter / 2) / increment) = mc.clust.parameter[p];
                    }
                }
            }
        }

        /* Write output */
        for (int p = 0; p < NPARAMS; p++)
        {
            if (ctrl.priorVar[p] > EPS || p == FEH || p == MOD || p == ABS)
            {
                ctrl.burninFile << boost::format("%10.6f ") % mc.clust.parameter[p];
            }
        }

        ctrl.burninFile << boost::format("%10.6f") % logPostCurr << endl;
    }

    ctrl.burninFile.close();

    make_cholesky_decomp(ctrl, params);

    // Main run
    ctrl.resFile.open(ctrl.clusterFilename);
    if (!ctrl.resFile)
    {
        cerr << "***Error: File " << ctrl.clusterFilename << " was not available for writing.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    printHeader (ctrl.resFile, ctrl.priorVar);

    for (int iteration = 0; iteration < ctrl.nIter * ctrl.thin; iteration++)
    {
        propClust = mc.clust;
        propClustCorrelated (propClust, ctrl);

        logPostProp = logPostStep (mc, evoModels, wdMass1Grid, propClust, fsLike, ltau);

        /* accept/reject */
        if (acceptClustMarg (logPostCurr, logPostProp, ltau))
        {
            mc.clust = propClust;
            logPostCurr = logPostProp;
            accept++;
        }
        else
        {
            reject++;
        }

        if (iteration % ctrl.thin == 0)
        {
            for (int p = 0; p < NPARAMS; p++)
            {
                if (ctrl.priorVar[p] > EPS || p == FEH || p == MOD || p == ABS)
                {
                    ctrl.resFile << boost::format("%10.6f ") % mc.clust.parameter[p];
                }
            }
            ctrl.resFile << boost::format("%10.6f") % logPostCurr << endl;
        }
    }

    ctrl.resFile.close();

    cout << "\nAcceptance ratio: " << (double) accept / (accept + reject) << endl;

    /* clean up */
    freeGlobalIso (isochrone);

    return 0;
}
