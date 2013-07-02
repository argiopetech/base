#include <array>
#include <atomic>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include <boost/format.hpp>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>

#include "densities.hpp"
#include "loadModels.hpp"
#include "marg.hpp"
#include "Model.hpp"
#include "mpiMcmc.hpp"
#include "mt19937ar.hpp"

using std::array;
using std::atomic;
using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::isfinite;

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
 * read control parameters from input stream
 */
void initIfmrMcmcControl (Chain &mc, struct ifmrMcmcControl &ctrl, const Model &evoModels, Settings &settings)
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

    /* set starting values for IFMR parameters */
    ctrl.priorMean[IFMR_SLOPE] = 0.08;
    ctrl.priorMean[IFMR_INTERCEPT] = 0.65;
    if (evoModels.IFMR <= 10)
        ctrl.priorMean[IFMR_QUADCOEF] = 0.0001;
    else
        ctrl.priorMean[IFMR_QUADCOEF] = 0.08;

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


/*
 * Initialize chain
 */
void initChain (Chain &mc, const struct ifmrMcmcControl &ctrl, const Model &evoModels, array<double, 2> &ltau)
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

        for (i = 0; i < ctrl.numFilts; i++)
        {
            star.photometry[i] = 0.0;
        }

        // find photometry for initial values of currentClust and mc.stars
        evolve (mc.clust, evoModels, star, ltau);

        if (star.status[0] == WD)
        {
            star.UStepSize = 0.05;      // use larger initial step size for white dwarfs
            star.massRatio = 0.0;
        }
    }
} /* initChain */


// Create Cholesky Decomp
void make_cholesky_decomp(struct ifmrMcmcControl &ctrl, Matrix<double, NPARAMS, nSave> &params)
{
    double cov;
    int nParamsUsed = 0;

    int nSave = 10;             /*changed from 100 to 10 */

    for (int p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar[p] > EPSILON)
        {
            nParamsUsed++;
        }
    }

    /* compute Cholesky decomposition of covariance matrix */
    int h, k;
    gsl_matrix *covMat = gsl_matrix_alloc (nParamsUsed, nParamsUsed);

    h = 0;

    double cholScale = 1000;    /* for numerical stability */

    for (int i = 0; i < NPARAMS; i++)
    {
        if (ctrl.priorVar[i] > EPSILON)
        {
            k = 0;
            for (int j = 0; j < NPARAMS; j++)
            {
                if (ctrl.priorVar[j] > EPSILON)
                {
                    cov = gsl_stats_covariance (params.at(i).data(), 1, params.at(j).data(), 1, nSave);
                    gsl_matrix_set (covMat, h, k, cov * cholScale * cholScale); /* for numerical stability? */

                    if (h != k)
                    {
                        gsl_matrix_set (covMat, k, h, cov * cholScale * cholScale);
                    }

                    k++;
                }
            }
            h++;
        }
    }

    cout << endl;

    for (int i = 0; i < nParamsUsed; i++)
    {
        for (int j = 0; j < nParamsUsed; j++)
        {
            cout << boost::format("%10.6f") % gsl_matrix_get (covMat, i, j) << " ";
        }
        cout << endl;
    }

    /* Cholesky decomposition */
    gsl_linalg_cholesky_decomp (covMat);

    /* compute proposal matrix from Cholesky factor */

    /* Gelman, Roberts, Gilks scale */
    double GRGscale = 0.97;     /* = 2.38 / sqrt(6) */

    h = 0;
    for (int i = 0; i < NPARAMS; i++)
    {
        if (ctrl.priorVar[i] > EPSILON)
        {
            k = 0;
            for (int j = 0; j < NPARAMS; j++)
            {
                if (ctrl.priorVar[j] > EPSILON)
                {
                    if (j <= i)
                    {
                        ctrl.propMatrix[i][j] = GRGscale * gsl_matrix_get (covMat, h, k) / cholScale;
                    }
                    else
                    {
                        ctrl.propMatrix[i][j] = 0.0;
                    }
                    k++;
                }
                else
                {
                    ctrl.propMatrix[i][j] = 0.0;
                }
            }
            h++;
        }
        else
        {
            for (int j = 0; j < NPARAMS; j++)
            {
                ctrl.propMatrix[i][j] = 0.0;
            }
        }
    }
}

void parallelFor(const unsigned int size, std::function<void(const unsigned int)> func)
{
    const unsigned int nbThreads = std::thread::hardware_concurrency();

    std::vector < std::thread > threads;

    for (unsigned int idThread = 0; idThread < nbThreads; idThread++) {
        auto threadFunc = [=, &threads]() {
            for (unsigned int i=idThread; i<size; i+=nbThreads) {
                func(i);
            }
        };

        threads.push_back(std::thread(threadFunc));
    }
    for (auto &t : threads)
        t.join();
}

double logPostStep(Chain &mc, const Model &evoModels, array<double, N_WD_MASS1> &wdMass1Grid, Cluster &propClust, double fsLike, array<double, 2> &ltau)
{
    atomic<double> postClusterStar(0.0);

    double logPostProp = logPriorClust (propClust, evoModels);

    if (isfinite(logPostProp))
    {
        /* loop over assigned stars */
        for (auto star : mc.stars)
        {
            /* loop over all (mass1, mass ratio) pairs */
            if (star.status[0] == WD)
            {
                postClusterStar = 0.0;

                for (int j = 0; j < N_WD_MASS1; j++)
                {
                    double tmpLogPost;
                    Star wd(star);
                    wd.boundsFlag = 0;
                    wd.isFieldStar = 0;
                    wd.U = wdMass1Grid[j];
                    wd.massRatio = 0.0;
                           
                    evolve (propClust, evoModels, wd, ltau);

                    if (wd.boundsFlag)
                    {
                        cerr <<"**wd[" << j << "].boundsFlag" << endl;
                    }
                    else
                    {
                        tmpLogPost = logPost1Star (wd, propClust, evoModels);
                        tmpLogPost += log ((mc.clust.M_wd_up - MIN_MASS1) / (double) N_WD_MASS1);

                        postClusterStar = postClusterStar +  exp (tmpLogPost);
                    }
                }
            }
            else
            {
                /* marginalize over isochrone */
                postClusterStar = margEvolveWithBinary (propClust, star, evoModels, ltau);
            }

            postClusterStar = postClusterStar * star.clustStarPriorDens;


            /* marginalize over field star status */
            logPostProp += log ((1.0 - star.clustStarPriorDens) * fsLike + postClusterStar);
        }
    }

    return logPostProp;
}

/* Decides whether to accept a proposed cluster property */
int acceptClustMarg (double logPostCurr, double logPostProp, array<double, 2> &ltau)
{
    if (isinf (logPostProp))
    {
        puts ("-Inf posterior proposed and rejected");
        return 0;
    }

    double alpha = logPostProp - logPostCurr;

    if (alpha >= 0)             // Short circuit exit to the MH algorithm
    {
        return 1;
    }

    double u = genrand_res53 ();

    if (u < 1.e-15)
        u = 1.e-15;
    u = log (u);

    if (u < alpha)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
