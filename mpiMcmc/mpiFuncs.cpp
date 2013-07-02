#include <array>
#include <atomic>
#include <iostream>
#include <thread>
#include <vector>

#include <boost/format.hpp>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>

#include "densities.hpp"
#include "marg.hpp"
#include "Model.hpp"
#include "mpiMcmc.hpp"
#include "mt19937ar.hpp"

using std::array;
using std::atomic;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::isfinite;

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
