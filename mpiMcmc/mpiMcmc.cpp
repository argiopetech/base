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
double ltau[2];
int aFilt = -1;

/* Used in densities.c. */
double filterPriorMin[FILTS];
double filterPriorMax[FILTS];
double priorMean[NPARAMS], priorVar[NPARAMS];
extern double ageLimit[2];      /* Defined in evolve.c, set in loadModels. */
extern struct globalIso isochrone;

int needMassNow = 0, useFilt[FILTS], numFilts = 0;

Settings settings;

/*
 * Initialize chain
 */
void initChain (Chain &mc, const struct ifmrMcmcControl &ctrl, Model &evoModels)
{
    int p;

    for (p = 0; p < NPARAMS; p++)
    {
        mc.acceptClust[p] = mc.rejectClust[p] = 0;
    }

    mc.clust.parameter[FEH] = ctrl.priorMean[FEH];
    mc.clust.parameter[MOD] = ctrl.priorMean[MOD];
    mc.clust.parameter[ABS] = ctrl.priorMean[ABS];
    mc.clust.parameter[YYY] = ctrl.priorMean[YYY];
    mc.clust.parameter[AGE] = ctrl.initialAge;
    mc.clust.mean[AGE] = ctrl.initialAge;
    mc.clust.mean[YYY] = ctrl.priorMean[YYY];
    mc.clust.mean[MOD] = ctrl.priorMean[MOD];
    mc.clust.mean[FEH] = ctrl.priorMean[FEH];
    mc.clust.mean[ABS] = ctrl.priorMean[ABS];
    mc.clust.betamabs = 0.0;
    mc.clust.betaFabs = 0.0;

    /* IFMR parameters */
    mc.clust.parameter[IFMR_SLOPE] = ctrl.priorMean[IFMR_SLOPE];
    mc.clust.parameter[IFMR_INTERCEPT] = ctrl.priorMean[IFMR_INTERCEPT];
    mc.clust.parameter[IFMR_QUADCOEF] = ctrl.priorMean[IFMR_QUADCOEF];
    mc.clust.mean[IFMR_SLOPE] = ctrl.priorMean[IFMR_SLOPE];
    mc.clust.mean[IFMR_INTERCEPT] = ctrl.priorMean[IFMR_INTERCEPT];
    mc.clust.mean[IFMR_QUADCOEF] = ctrl.priorMean[IFMR_QUADCOEF];


    int i;

    for (auto star : mc.stars)
    {
        star.meanMassRatio = 0.0;
        star.isFieldStar = 0;
        star.clustStarProposalDens = star.clustStarPriorDens;   // Use prior prob of being clus star
        star.UStepSize = 0.001; // within factor of ~2 for most main sequence stars
        star.massRatioStepSize = 0.001;

        for (i = 0; i < NPARAMS; i++)
        {
            star.beta[i][0] = 0.0;
            star.beta[i][1] = 0.0;
        }

        star.betaMassRatio[0] = 0.0;
        star.betaMassRatio[1] = 0.0;
        star.meanU = 0.0;
        star.varU = 0.0;

        for (i = 0; i < 2; i++)
            star.wdType[i] = 0;

        for (i = 0; i < numFilts; i++)
        {
            star.photometry[i] = 0.0;
        }

        // find photometry for initial values of currentClust and mc.stars
        evolve (&mc.clust, evoModels, star);

        if (star.status[0] == WD)
        {
            star.UStepSize = 0.05;      // use larger initial step size for white dwarfs
            star.massRatio = 0.0;
        }
    }
} /* initChain */



/*
 * read control parameters from input stream
 */
void initIfmrMcmcControl (Chain &mc, struct ifmrMcmcControl &ctrl, Model &evoModels)
{

    double priorSigma;

    ctrl.numFilts = 0;

    int ii;

    for (ii = 0; ii < FILTS; ii++)
        ctrl.useFilt[ii] = 0;

    /* Read number of steps, burn-in details, random seed */
    init_genrand (settings.seed);

    /* load models */
    loadModels (&mc.clust, evoModels, settings);

    ctrl.priorMean[FEH] = settings.cluster.Fe_H;
    priorSigma = settings.cluster.sigma.Fe_H;

    if (priorSigma < 0.0)
    {
        priorSigma = 0.0;
    }
    ctrl.priorVar[FEH] = priorSigma * priorSigma;

    ctrl.priorMean[MOD] = settings.cluster.distMod;
    priorSigma = settings.cluster.sigma.distMod;

    if (priorSigma < 0.0)
    {
        priorSigma = 0.0;
    }
    ctrl.priorVar[MOD] = priorSigma * priorSigma;

    ctrl.priorMean[ABS] = settings.cluster.Av;
    priorSigma = settings.cluster.sigma.Av;

    if (priorSigma < 0.0)
    {
        priorSigma = 0.0;
    }
    ctrl.priorVar[ABS] = priorSigma * priorSigma;

    ctrl.initialAge = settings.cluster.logClusAge;
    ctrl.priorVar[AGE] = 1.0;

    if (evoModels.IFMR <= 3)
    {
        ctrl.priorVar[IFMR_SLOPE] = 0.0;
        ctrl.priorVar[IFMR_INTERCEPT] = 0.0;
        ctrl.priorVar[IFMR_QUADCOEF] = 0.0;
    }
    else if (evoModels.IFMR <= 8)
    {
        ctrl.priorVar[IFMR_SLOPE] = 1.0;
        ctrl.priorVar[IFMR_INTERCEPT] = 1.0;
        ctrl.priorVar[IFMR_QUADCOEF] = 0.0;
    }
    else
    {
        ctrl.priorVar[IFMR_SLOPE] = 1.0;
        ctrl.priorVar[IFMR_INTERCEPT] = 1.0;
        ctrl.priorVar[IFMR_QUADCOEF] = 1.0;
    }

    // copy values to global variables
    priorVar[AGE] = ctrl.priorVar[AGE];
    priorVar[FEH] = ctrl.priorVar[FEH];
    priorVar[MOD] = ctrl.priorVar[MOD];
    priorVar[ABS] = ctrl.priorVar[ABS];
    priorVar[IFMR_SLOPE] = ctrl.priorVar[IFMR_SLOPE];
    priorVar[IFMR_INTERCEPT] = ctrl.priorVar[IFMR_INTERCEPT];
    priorVar[IFMR_QUADCOEF] = ctrl.priorVar[IFMR_QUADCOEF];

    priorMean[FEH] = ctrl.priorMean[FEH];
    priorMean[MOD] = ctrl.priorMean[MOD];
    priorMean[ABS] = ctrl.priorMean[ABS];

    /* set starting values for IFMR parameters */
    ctrl.priorMean[IFMR_SLOPE] = 0.08;
    ctrl.priorMean[IFMR_INTERCEPT] = 0.65;
    if (evoModels.IFMR <= 10)
        ctrl.priorMean[IFMR_QUADCOEF] = 0.0001;
    else
        ctrl.priorMean[IFMR_QUADCOEF] = 0.08;
    priorMean[IFMR_SLOPE] = ctrl.priorMean[IFMR_SLOPE];
    priorMean[IFMR_INTERCEPT] = ctrl.priorMean[IFMR_INTERCEPT];
    priorMean[IFMR_QUADCOEF] = ctrl.priorMean[IFMR_QUADCOEF];

    /* open model file, choose model set, and load models */
    if (settings.mainSequence.msRgbModel == MsModel::CHABHELIUM)
    {
        scanf ("%lf %lf", &ctrl.priorMean[YYY], &ctrl.priorVar[YYY]);

        if (ctrl.priorVar[YYY] < 0.0)
        {
            ctrl.priorVar[YYY] = 0.0;
        }
    }
    else
    {
        ctrl.priorMean[YYY] = 0.0;
        ctrl.priorVar[YYY] = 0.0;
    }
    priorVar[YYY] = ctrl.priorVar[YYY];
    priorMean[YYY] = ctrl.priorMean[YYY];

    /* read burnIter and nIter */
    ctrl.burnIter = settings.mpiMcmc.burnIter;
    ctrl.nIter = settings.mpiMcmc.maxIter;
    ctrl.thin = settings.mpiMcmc.thin;

    /* open files for reading (data) and writing */
    string filename;

    filename = settings.files.phot;
    ctrl.rData.open(filename);
    if (!ctrl.rData)
    {
        cerr << "***Error: Photometry file " << filename << " was not found.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    ctrl.minMag = settings.cluster.minMag;
    ctrl.maxMag = settings.cluster.maxMag;
    ctrl.iMag = settings.cluster.index;

    if (ctrl.iMag < 0 || ctrl.iMag > FILTS)
    {
        cerr << "***Error: " << ctrl.iMag << " not a valid magnitude index.  Choose 0, 1,or 2.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    ctrl.clusterFilename = settings.files.output + ".res";

    ctrl.iStart = 0;

    /* Initialize filter prior mins and maxes */
    int j;

    for (j = 0; j < FILTS; j++)
    {
        ctrl.filterPriorMin[j] = 1000;
        ctrl.filterPriorMax[j] = -1000;
    }
} /* initIfmrMcmcControl */



void initMassGrids (array<double, N_MS_MASS1 * N_MS_MASS_RATIO> &msMass1Grid, array<double, N_MS_MASS1 * N_MS_MASS_RATIO> &msMassRatioGrid, array<double, N_WD_MASS1> &wdMass1Grid, Chain const &mc)
{
    double maxMass1 = mc.clust.M_wd_up;
    double mass1, massRatio;
    double dMsMass1 = (maxMass1 - MIN_MASS1) / (double) N_MS_MASS1;
    double dMsMassRatio = 1.0 / (double) N_MS_MASS_RATIO;
    double dWdMass1 = (maxMass1 - MIN_MASS1) / (double) N_WD_MASS1;

    int i = 0;

    for (mass1 = MIN_MASS1; mass1 < maxMass1; mass1 += dMsMass1)
    {
        for (massRatio = 0.0; massRatio < 1.0; massRatio += dMsMassRatio)
        {
            msMass1Grid[i] = mass1;
            msMassRatioGrid[i] = massRatio;
            i++;
        }
    }

    i = 0;
    for (mass1 = MIN_MASS1; mass1 < maxMass1; mass1 += dWdMass1)
    {
        wdMass1Grid[i] = mass1;
        i++;
    }
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

    array<double, N_MS_MASS1 * N_MS_MASS_RATIO> msMass1Grid;
    array<double, N_MS_MASS1 * N_MS_MASS_RATIO> msMassRatioGrid;
    array<double, N_WD_MASS1> wdMass1Grid;

    Matrix<double, NPARAMS, nSave> params;

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

    Model evoModels = makeModel(settings);

    increment = settings.mpiMcmc.burnIter / (2 * nSave);

    initIfmrMcmcControl (mc, ctrl, evoModels);

    evoModels.WDatm = BERGERON;

    for (int p = 0; p < NPARAMS; p++)
    {
        mc.clust.priorVar[p] = ctrl.priorVar[p];
        mc.clust.priorMean[p] = ctrl.priorMean[p];
    }

    readCmdData (mc, ctrl, evoModels);

    evoModels.numFilts = ctrl.numFilts;
    numFilts = ctrl.numFilts;

    initChain (mc, ctrl, evoModels);

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

        logPostProp = logPostStep (mc, evoModels, wdMass1Grid, propClust, fsLike);

        /* accept/reject */
        if (acceptClustMarg (logPostCurr, logPostProp))
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

        logPostProp = logPostStep (mc, evoModels, wdMass1Grid, propClust, fsLike);

        /* accept/reject */
        if (acceptClustMarg (logPostCurr, logPostProp))
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
    freeGlobalIso (&isochrone);

    return 0;
}
