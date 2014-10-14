#ifndef CHAIN_HPP
#define CHAIN_HPP

#include <array>
#include <ostream>
#include <iostream>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>

#include "Cluster.hpp"
#include "Matrix.hpp"
#include "McmcApplication.hpp"
#include "Utility.hpp"

using std::ostream;
using std::cout;
using std::endl;
using std::cerr;

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
    {}

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
        int nParamsUsed = 0;

        Matrix<double, 2, 2> propMatrix;

        /* compute Cholesky decomposition of covariance matrix */
        gsl_matrix *covMat = gsl_matrix_alloc (2, 2);

        double cholScale = 1000;    /* for numerical stability */

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
        gsl_linalg_cholesky_decomp (covMat);

        /* compute proposal matrix from Cholesky factor */

        /* Gelman, Roberts, Gilks scale */
        double GRGscale = 2.38 / sqrt(nParamsUsed);

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
