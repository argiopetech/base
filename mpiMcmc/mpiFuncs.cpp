#include <array>
#include <iostream>
#include <vector>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>

#include "densities.hpp"
#include "marg.hpp"
#include "mpiMcmc.hpp"
#include "mt19937ar.hpp"

using std::array;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

// Create Cholesky Decomp
void make_cholesky_decomp(struct ifmrMcmcControl &ctrl, double **params)
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
                    cov = gsl_stats_covariance (params[i], 1, params[j], 1, nSave);
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

    for (int i = 0; i < nParamsUsed; i++)
    {
        for (int j = 0; j < nParamsUsed; j++)
        {
            cout << gsl_matrix_get (covMat, i, j) << " ";
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

double logPostStep(struct chain &mc, array<double, N_WD_MASS1> &wdMass1Grid, struct cluster &propClust, double fsLike)
{
    vector<struct star> wd(N_WD_MASS1);
    double postClusterStar = 0.0;

    double logPostProp = logPriorClust (&propClust);

    if (isfinite(logPostProp))
    {
        /* loop over assigned stars */
        for (int i = 0; i < mc.clust.nStars; i++)
        {
            /* loop over all (mass1, mass ratio) pairs */
            if (mc.stars[i].status[0] == WD)
            {

                postClusterStar = 0.0;
                double tmpLogPost;

                for (int j = 0; j < N_WD_MASS1; j++)
                {
                    wd[j] = mc.stars[i];
                    wd[j].boundsFlag = 0;
                    wd[j].isFieldStar = 0;
                    wd[j].U = wdMass1Grid[j];
                    wd[j].massRatio = 0.0;

                    evolve (&propClust, wd, j);

                    if (wd[j].boundsFlag)
                    {
                        cerr <<"**wd[" << j << "].boundsFlag" << endl;
                    }
                    else
                    {
                        tmpLogPost = logPost1Star (&wd[j], &propClust);
                        tmpLogPost += log ((mc.clust.M_wd_up - MIN_MASS1) / (double) N_WD_MASS1);

                        postClusterStar += exp (tmpLogPost);
                    }
                }
            }
            else
            {
                /* marginalize over isochrone */
                postClusterStar = margEvolveWithBinary (&propClust, &mc.stars[i]);
            }

            postClusterStar *= mc.stars[i].clustStarPriorDens;


            /* marginalize over field star status */
            logPostProp += log ((1.0 - mc.stars[i].clustStarPriorDens) * fsLike + postClusterStar);
        }
    }

    return logPostProp;
}

/* Decides whether to accept a proposed cluster property */
int acceptClustMarg (double logPostCurr, double logPostProp)
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
