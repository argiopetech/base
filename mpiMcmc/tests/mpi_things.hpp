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
#include "../mpiMcmc.hpp"
#include "../mpiFuncs.hpp"
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
int aFilt = -1;

/* Used in densities.c. */
double filterPriorMin[FILTS];
double filterPriorMax[FILTS];
extern double ageLimit[2];      /* Defined in evolve.c, set in loadModels. */
extern struct globalIso isochrone;

int needMassNow = 0, useFilt[FILTS], numFilts = 0;

Settings settings;


void initStepSizes (Cluster &clust)
{
    clust.stepSize[AGE] = 0.005;
    clust.stepSize[FEH] = 0.005;
    clust.stepSize[MOD] = 0.005;
    clust.stepSize[ABS] = 0.002;
    clust.stepSize[YYY] = 0.002;
    clust.stepSize[IFMR_INTERCEPT] = 0.01;
    clust.stepSize[IFMR_SLOPE] = 0.008;
    clust.stepSize[IFMR_QUADCOEF] = 0.008;
}



/*
 * Read data
 */
void readCmdData (Chain &mc, struct ifmrMcmcControl &ctrl, Model &evoModels)
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
                evoModels.numFilts++;
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
    numFilts = ctrl.numFilts;

    assert(mc.clust.nStars == mc.stars.size());
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
