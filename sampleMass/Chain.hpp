#ifndef CHAIN_HPP
#define CHAIN_HPP

#include <array>
#include <ostream>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>

#include "Cluster.hpp"
#include "Matrix.hpp"
#include "McmcApplication.hpp"
#include "StellarModel.hpp"
#include "Utility.hpp"

using std::ostream;
using std::cout;
using std::endl;
using std::cerr;

class NonPositiveDefiniteMatrix : public std::domain_error
{
  public:
    explicit NonPositiveDefiniteMatrix (const std::string& what_arg)
        : std::domain_error(what_arg) {}

    explicit NonPositiveDefiniteMatrix (const char* what_arg)
        : std::domain_error(what_arg) {}
};

template <class T>
class Chain : public McmcApplication
{
  private:
    MVatrix<double, 2> params;

    T clust;

    ostream &fout;

    double logPostCurr = -std::numeric_limits<double>::infinity();

  public:
    Chain(uint32_t seed, T clust, ostream &fout)
        : McmcApplication(seed), clust(clust), fout(fout)
    {
        // Disable the GSL error handler. This should be able to be
        // called multiple times without adverse effect. After this,
        // all GSL functions will return success/failure in an integer
        // rather than exiting.
        gsl_set_error_handler_off();
    }

    void reset()
    {
        resetRatio();

        for ( auto &p : params )
        {
            p.clear();
        }
    }

    void run(std::function<T(T)> propose, std::function<double(T&)> logPost, std::function<void(const T&)> checkPriors, int iters, int thin = 1)
    {
        for (int iteration = 0; iteration < iters * thin; iteration++)
        {
            double logPostProp = 0.0;
            auto propClust = propose (clust);

            try
            {
                checkPriors(propClust);
                logPostProp = logPost (propClust);
            }
            catch(InvalidCluster &e)
            {
                logPostProp = -std::numeric_limits<double>::infinity();
            }
            catch(InvalidModelError &e)
            {
                logPostProp = -std::numeric_limits<double>::infinity();
            }

            /* accept/reject */
            if (acceptClustMarg (logPostCurr, logPostProp))
            {
                clust = propClust;
                logPostCurr = logPostProp;
            }

            params.at(0).push_back(clust.primary.mass);
            params.at(1).push_back(clust.secondary.mass);
        }
    }

    Matrix<double, 2, 2> makeCholeskyDecomp() const
    {
        double cov;

        Matrix<double, 2, 2> propMatrix;

        /* compute Cholesky decomposition of covariance matrix */
        gsl_matrix *covMat = gsl_matrix_alloc (2, 2);

        const double cholScale = 1000;    /* for numerical stability */

        {
            cov = gsl_stats_covariance (params.at(0).data(), 1, params.at(0).data(), 1, params.at(0).size());
            gsl_matrix_set (covMat, 0, 0, cov * cholScale * cholScale);
        }

        {
            cov = gsl_stats_covariance (params.at(1).data(), 1, params.at(0).data(), 1, params.at(1).size());
            gsl_matrix_set (covMat, 1, 0, cov * cholScale * cholScale);
            gsl_matrix_set (covMat, 0, 1, cov * cholScale * cholScale);
        }

        {
            cov = gsl_stats_covariance (params.at(1).data(), 1, params.at(1).data(), 1, params.at(1).size());
            gsl_matrix_set (covMat, 1, 1, cov * cholScale * cholScale);
        }

        /* Cholesky decomposition */
        if (gsl_linalg_cholesky_decomp (covMat))
        {
            // Handle a potential error.
            // As far as I know, this will always be due to a
            // non-positive definite matrix.
            throw NonPositiveDefiniteMatrix("Non-positive matrix in makeCholeskyDecomp");
        }

        /* compute proposal matrix from Cholesky factor */

        /* Gelman, Roberts, Gilks scale */
        double GRGscale = 2.38 / sqrt(2);

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                if (j <= i)
                {
                    propMatrix.at(i).at(j) = GRGscale * gsl_matrix_get (covMat, i, j) / cholScale;
                }
            }
        }

        gsl_matrix_free (covMat);

        return propMatrix;
    }

    T get() const { return clust; };
};

#endif
