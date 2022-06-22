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
#include "IO/BackingStore.hpp"
#include "IO/Records.hpp"
#include "Matrix.hpp"
#include "McmcApplication.hpp"

using std::cout;
using std::endl;
using std::cerr;

template <class T>
class Chain : public McmcApplication
{
  private:
    MVatrix<double, NPARAMS> params;
    std::array<double, NPARAMS> priorVar;

    T clust;

    SinglePopBackingStore &store;
    const bool modIsParallax;

    double logPostCurr = -std::numeric_limits<double>::infinity();

  public:
    Chain(uint32_t seed, std::array<double, NPARAMS> priorVar, T clust, SinglePopBackingStore &store, bool modIsParallax)
        : McmcApplication(seed), priorVar(priorVar), clust(clust), store(store), modIsParallax(modIsParallax)
    {}

    void reset()
    {
        resetRatio();

        for ( auto &p : params )
        {
            p.clear();
        }
    }

    void run(AdaptiveMcmcStage stage, std::function<T(T)> propose, std::function<double(T&)> logPost, std::function<void(const T&)> checkPriors, int iters, int thin = 1)
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
            catch(MsBoundsError &e)
            {
                logPostProp = -std::numeric_limits<double>::infinity();
            }

            /* accept/reject */
            if (acceptClustMarg (logPostCurr, logPostProp))
            {
                clust = propClust;
                logPostCurr = logPostProp;
            }

            /* save draws to estimate covariance matrix for more efficient Metropolis */
            for (int p = 0; p < NPARAMS; p++)
            {
                if (priorVar[p] > EPSILON)
                {
                    params.at(p).push_back(clust.getParam(p));
                }
            }

            if (iteration % thin == 0)
            {
                store.save({store.nextIteration(), stage, clust, logPostCurr, modIsParallax});
            }
        }
    }

    Matrix<double, NPARAMS, NPARAMS> makeCholeskyDecomp() const
    {
        double cov;
        int nParamsUsed = 0;

        Matrix<double, NPARAMS, NPARAMS> propMatrix;

        for (int p = 0; p < NPARAMS; p++)
        {
            if (priorVar.at(p) > EPSILON)
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
            if (priorVar.at(i) > EPSILON)
            {
                k = 0;
                for (int j = 0; j < NPARAMS; j++)
                {
                    if (priorVar.at(j) > EPSILON)
                    {
                        cov = gsl_stats_covariance (params.at(i).data(), 1, params.at(j).data(), 1, params.at(i).size());
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

        // for (int k = 0; k < params.at(0).size(); ++k)
        // {
        //     for (int i = 0; i < NPARAMS; i++)
        //     {
        //         if (priorVar.at(i) > EPSILON || i == YYA || i == YYB || i == LAMBDA)
        //         {
        //             cout << base::utility::format("%10.6f") % params.at(i).at(k) << " ";
        //         }
        //     }

        //     cout << '\n';
        // }

        // cout << endl;

        // for (int i = 0; i < nParamsUsed; i++)
        // {
        //     for (int j = 0; j < nParamsUsed; j++)
        //     {
        //         cout << base::utility::format("%10.6f") % gsl_matrix_get (covMat, i, j) << " ";
        //     }
        //     cout << endl;
        // }

        /* Cholesky decomposition */
        gsl_linalg_cholesky_decomp (covMat);

        /* compute proposal matrix from Cholesky factor */

        /* Gelman, Roberts, Gilks scale */
        double GRGscale = 2.38 / sqrt(nParamsUsed);

        h = 0;
        for (int i = 0; i < NPARAMS; i++)
        {
            if (priorVar.at(i) > EPSILON)
            {
                k = 0;
                for (int j = 0; j < NPARAMS; j++)
                {
                    if (priorVar.at(j) > EPSILON)
                    {
                        if (j <= i)
                        {
                            propMatrix.at(i).at(j) = GRGscale * gsl_matrix_get (covMat, h, k) / cholScale;
                        }
                        else
                        {
                            propMatrix.at(i).at(j) = 0.0;
                        }
                        k++;
                    }
                    else
                    {
                        propMatrix.at(i).at(j) = 0.0;
                    }
                }
                h++;
            }
            else
            {
                for (int j = 0; j < NPARAMS; j++)
                {
                    propMatrix.at(i).at(j) = 0.0;
                }
            }
        }

        return propMatrix;
    }

    T getCluster() const { return clust; };
};

#endif
