#include <array>
#include <iostream>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>

#include "mpiMcmc.hpp"

using std::array;
using std::cout;
using std::endl;

// Create Cholesky Decomp
void make_cholesky_decomp(struct ifmrMcmcControl &ctrl, double **params)
{
    double cov;
    int nParamsUsed = 0;

    int nSave = 10;             /*changed from 100 to 10 */
    int increment = ctrl.burnIter / (2 * nSave);

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
