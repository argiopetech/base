#ifndef CHAIN_HPP
#define CHAIN_HPP

#include <array>
#include <fstream>
#include <iostream>

#include "Cluster.hpp"
#include "Matrix.hpp"
#include "McmcApplication.hpp"

using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

class Chain : public McmcApplication
{
  private:
    MVatrix<double, NPARAMS> params;
    std::array<double, NPARAMS> priorVar;

    Cluster clust;

    ofstream &fout;

    double logPostCurr = std::numeric_limits<double>::lowest();

  public:
    Chain(uint32_t seed, std::array<double, NPARAMS> priorVar, Cluster clust, ofstream &fout)
        : McmcApplication(seed), priorVar(priorVar), clust(clust), fout(fout)
    {}

    void reset()
    {
        resetRatio();

        for ( auto &p : params )
        {
            p.clear();
        }
    }

    void run(std::function<Cluster(Cluster)> propose, std::function<double(Cluster&)> logPost, int iters, int thin = 1)
    {
        for (int iteration = 0; iteration < iters * thin; iteration++)
        {
            double logPostProp = 0.0;
            auto propClust = propose (clust);

            try
            {
                logPostProp = logPost (propClust);
            }
            catch(InvalidCluster &e)
            {
                logPostProp = std::numeric_limits<double>::infinity();
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
                for (int p = 0; p < NPARAMS; p++)
                {
                    if (priorVar.at(p) > EPS || p == FEH || p == MOD || p == ABS)
                    {
                        fout << boost::format("%10.6f ") % clust.getParam(p);
                    }
                }
                fout << boost::format("%10.6f") % logPostCurr << endl;
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

    Cluster getCluster() const { return clust; };
};

#endif
