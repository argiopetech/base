#include <array>
#include <iostream>
#include <fstream>
#include <mutex>
#include <string>
#include <sstream>
#include <thread>
#include <vector>
#include <stdexcept>

#include <boost/format.hpp>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>

#include "densities.hpp"
#include "MsFilterSet.hpp"
#include "loadModels.hpp"
#include "marg.hpp"
#include "Model.hpp"
#include "mpiMcmc.hpp"
#include "mt19937ar.hpp"
#include "samplers.hpp"
#include "Utility.hpp"
#include "WhiteDwarf.hpp"

using std::array;
using std::mutex;
using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::isfinite;
using std::ofstream;
using std::istringstream;

using base::utility::ThreadPool;

/*
 * Read data
 */
void readCmdData (vector<Star> &stars, struct ifmrMcmcControl &ctrl, const Model &evoModels, vector<int> &filters, std::array<double, FILTS> &filterPriorMin, std::array<double, FILTS> &filterPriorMax, const Settings &settings)
{
    string line, pch;

    //Parse the header of the file to determine which filters are being used
    getline(ctrl.rData, line);  // Read in the header line

    istringstream header(line); // Ignore the first token (which is "id")

    header >> pch;

    while (!header.eof())
    {
        header >> pch;

        if (pch == "sig")
            break;                      // and check to see if they are 'sig'.  If they are, there are no more filters

        for (int filt = 0; filt < FILTS; filt++)
        {                               // Otherwise check to see what this filter's name is
            if (pch == evoModels.filterSet->getFilterName(filt))
            {
                filters.push_back(filt);
                ctrl.numFilts++;
                const_cast<Model&>(evoModels).numFilts++;
                break;
            }
        }
    }

    // This loop reads in photometry data
    // It also reads a best guess for the mass
    stars.clear();

    while (!ctrl.rData.eof())
    {
        getline(ctrl.rData, line);

        if (ctrl.rData.eof())
            break;

        stars.push_back(Star(line, ctrl.numFilts));

        for (int i = 0; i < ctrl.numFilts; i++)
        {
            if (stars.back().obsPhot.at(i) < filterPriorMin.at(i))
            {
                filterPriorMin.at(i) = stars.back().obsPhot.at(i);
            }

            if (stars.back().obsPhot.at(i) > filterPriorMax.at(i))
            {
                filterPriorMax.at(i) = stars.back().obsPhot.at(i);
            }
        }

        if (!(stars.back().status.at(0) == 3 || (stars.back().obsPhot.at(settings.cluster.index) >= settings.cluster.minMag && stars.back().obsPhot.at(settings.cluster.index) <= settings.cluster.maxMag)))
        {
            stars.pop_back();
        }
    }
} /* readCmdData */

void propClustBigSteps (Cluster &clust, struct ifmrMcmcControl const &ctrl)
{
    /* DOF defined in densities.h */
    double scale = 5.0;
    int p;

    for (p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar.at(p) > EPSILON)
        {
            clust.parameter.at(p) += sampleT (scale * scale * clust.stepSize.at(p) * clust.stepSize.at(p));
        }
    }
}

void propClustIndep (Cluster &clust, struct ifmrMcmcControl const &ctrl)
{
    /* DOF defined in densities.h */
    int p;

    for (p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar.at(p) > EPSILON)
        {
            clust.parameter.at(p) += sampleT (clust.stepSize.at(p) * clust.stepSize.at(p));
        }
    }
}

void propClustCorrelated (Cluster &clust, struct ifmrMcmcControl const &ctrl)
{
    /* DOF defined in densities.h */
    array<double, NPARAMS> indepProps;
    array<double, NPARAMS> corrProps;

    indepProps.fill(0.0);
    corrProps.fill(0.0);

    int p, k;

    for (p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar.at(p) > EPSILON)
        {
            indepProps.at(p) = sampleT (1.0);
        }
    }
    for (p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar.at(p) > EPSILON)
        {
            for (k = 0; k <= p; k++)
            {                           /* propMatrix is lower diagonal */
                if (ctrl.priorVar.at(k) > EPSILON)
                {
                    corrProps.at(p) += ctrl.propMatrix.at(p).at(k) * indepProps.at(k);
                }
            }
            clust.parameter.at(p) += corrProps.at(p);
        }
    }
}


void printHeader (ofstream &file, array<double, NPARAMS> const &priors)
{
    const array<string, NPARAMS> paramNames = { "    logAge",
                                                "         Y",
                                                "       FeH",
                                                "   modulus",
                                                "absorption",
                                                " IFMRconst",
                                                "   IFMRlin",
                                                "  IFMRquad"};

    for (int p = 0; p < NPARAMS; p++)
    {
        if (priors.at(p) > EPSILON || p == MOD || p == FEH || p == ABS)
        {
            file << paramNames.at(p) << ' ';
        }
    }
    file << "logPost" << endl;
}

void initMassGrids (array<double, N_MS_MASS1 * N_MS_MASS_RATIO> &msMass1Grid, array<double, N_MS_MASS1 * N_MS_MASS_RATIO> &msMassRatioGrid, array<double, N_WD_MASS1> &wdMass1Grid, Chain const &mc)
{
    double dMsMass1 = (mc.clust.M_wd_up - MIN_MASS1) / (double) N_MS_MASS1;
    double dMsMassRatio = 1.0 / (double) N_MS_MASS_RATIO;
    double dWdMass1 = (mc.clust.M_wd_up - MIN_MASS1) / (double) N_WD_MASS1;

    for (int i = 0; i < N_MS_MASS1; ++i)
    {
        for (int j = 0; j < N_MS_MASS_RATIO; ++j)
        {
            msMass1Grid.at(i) = MIN_MASS1 + (dMsMass1 * i);
            msMassRatioGrid.at(i) = dMsMassRatio * j;
        }
    }

    for (int i = 0; i < N_WD_MASS1; ++i)
    {
        wdMass1Grid.at(i) = MIN_MASS1 + (dWdMass1 * i);
    }
}

/*
 * read control parameters from input stream
 */
void initIfmrMcmcControl (Cluster &clust, struct ifmrMcmcControl &ctrl, const Model &evoModels, Settings &settings)
{
    ctrl.numFilts = 0;

    /* Read number of steps, burn-in details, random seed */
    init_genrand (settings.seed);

    ctrl.priorVar.fill(0);
    clust.parameter.fill(0);

    /* load models */
    loadModels (clust, evoModels, settings);

    clust.parameter[FEH] = settings.cluster.Fe_H;
    ctrl.priorVar[FEH] = settings.cluster.sigma.Fe_H;

    clust.parameter[MOD] = settings.cluster.distMod;
    ctrl.priorVar[MOD] = settings.cluster.sigma.distMod;

    clust.parameter[ABS] = settings.cluster.Av;
    ctrl.priorVar[ABS] = settings.cluster.sigma.Av;

    clust.parameter[AGE] = settings.cluster.logClusAge;
    ctrl.priorVar[AGE] = 1.0;

    if (settings.mainSequence.msRgbModel == MsModel::CHABHELIUM)
    {
        clust.parameter[YYY] = settings.cluster.Y;
        ctrl.priorVar[YYY] = settings.cluster.sigma.Y;
    }
    else
    {
        clust.parameter[YYY] = 0.0;
        ctrl.priorVar[YYY] = 0.0;
    }


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
    clust.parameter[IFMR_SLOPE] = 0.08;
    clust.parameter[IFMR_INTERCEPT] = 0.65;

    if (evoModels.IFMR <= 10)
        clust.parameter[IFMR_QUADCOEF] = 0.0001;
    else
        clust.parameter[IFMR_QUADCOEF] = 0.08;

    for (auto &var : ctrl.priorVar)
    {
        if (var < 0.0)
        {
            var = 0.0;
        }
        else
        {
            var = var * var;
        }
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

    if (settings.cluster.index < 0 || settings.cluster.index > FILTS)
    {
        cerr << "***Error: " << settings.cluster.index << " not a valid magnitude index.  Choose 0, 1,or 2.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    ctrl.clusterFilename = settings.files.output + ".res";

    std::copy(clust.parameter.begin(), clust.parameter.end(), clust.mean.begin());
    std::copy(clust.parameter.begin(), clust.parameter.end(), clust.priorMean.begin());
    std::copy(ctrl.priorVar.begin(), ctrl.priorVar.end(), clust.priorVar.begin());
} /* initIfmrMcmcControl */


/*
 * Initialize chain
 */
void initChain (Chain &mc, const Model &evoModels, array<double, 2> &ltau, const vector<int> &filters)
{
    array<double, FILTS> globalMags;

    for (int p = 0; p < NPARAMS; p++)
    {
        mc.acceptClust.at(p) = mc.rejectClust.at(p) = 0;
    }

    for (auto star : mc.stars)
    {
        star.clustStarProposalDens = star.clustStarPriorDens;   // Use prior prob of being clus star
        star.UStepSize = 0.001; // within factor of ~2 for most main sequence stars
        star.massRatioStepSize = 0.001;

        // find photometry for initial values of currentClust and mc.stars
        mc.clust.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, mc.clust.getFeH(), mc.clust.getY(), mc.clust.getAge());    // determine AGBt ZAMS mass, to find evol state
        evolve (mc.clust, evoModels, globalMags, filters, star, ltau);

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
        if (ctrl.priorVar.at(p) > EPSILON)
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
        if (ctrl.priorVar.at(i) > EPSILON)
        {
            k = 0;
            for (int j = 0; j < NPARAMS; j++)
            {
                if (ctrl.priorVar.at(j) > EPSILON)
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
        if (ctrl.priorVar.at(i) > EPSILON)
        {
            k = 0;
            for (int j = 0; j < NPARAMS; j++)
            {
                if (ctrl.priorVar.at(j) > EPSILON)
                {
                    if (j <= i)
                    {
                        ctrl.propMatrix.at(i).at(j) = GRGscale * gsl_matrix_get (covMat, h, k) / cholScale;
                    }
                    else
                    {
                        ctrl.propMatrix.at(i).at(j) = 0.0;
                    }
                    k++;
                }
                else
                {
                    ctrl.propMatrix.at(i).at(j) = 0.0;
                }
            }
            h++;
        }
        else
        {
            for (int j = 0; j < NPARAMS; j++)
            {
                ctrl.propMatrix.at(i).at(j) = 0.0;
            }
        }
    }
}

static ThreadPool pool(2);

double logPostStep(Chain &mc, const Model &evoModels, array<double, N_WD_MASS1> &wdMass1Grid, Cluster &propClust, double fsLike, const vector<int> &filters, std::array<double, FILTS> &filterPriorMin, std::array<double, FILTS> &filterPriorMax)
{
    mutex logPostMutex;

    double logPostProp = logPriorClust (propClust, evoModels);

    auto stars = mc.stars;

    propClust.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, propClust.getFeH(), propClust.getY(), propClust.getAge());    // determine AGBt ZAMS mass, to find evol state

    if (isfinite(logPostProp))
    {
        /* loop over assigned stars */
        pool.parallelFor(mc.stars.size(), [=,&logPostMutex,&logPostProp](int i)
        {
            double postClusterStar = 0.0;

            array<double, FILTS> globalMags;
            array<double, 2> ltau;

            /* loop over all (mass1, mass ratio) pairs */
            if (stars.at(i).status[0] == WD)
            {
                postClusterStar = 0.0;

                for (int j = 0; j < N_WD_MASS1; j++)
                {
                    double tmpLogPost;
                    Star wd(stars.at(i));
                    wd.isFieldStar = 0;
                    wd.U = wdMass1Grid.at(j);
                    wd.massRatio = 0.0;

                    evolve (propClust, evoModels, globalMags, filters, wd, ltau);

                    try
                    {
                        tmpLogPost = logPost1Star (wd, propClust, evoModels, filterPriorMin, filterPriorMax);
                        tmpLogPost += log ((mc.clust.M_wd_up - MIN_MASS1) / (double) N_WD_MASS1);

                        postClusterStar = postClusterStar +  exp (tmpLogPost);
                    }
                    catch ( WDBoundsError &e )
                    {
                        cerr << e.what() << endl;
                    }
                }
            }
            else
            {
                /* marginalize over isochrone */
                postClusterStar = margEvolveWithBinary (propClust, stars.at(i), evoModels, filters, ltau, globalMags, filterPriorMin, filterPriorMax);
            }

            postClusterStar = postClusterStar * stars.at(i).clustStarPriorDens;


            /* marginalize over field star status */
            std::lock_guard<mutex> lk(logPostMutex);
            logPostProp += log ((1.0 - stars.at(i).clustStarPriorDens) * fsLike + postClusterStar);
        });
    }

    return logPostProp;
}
