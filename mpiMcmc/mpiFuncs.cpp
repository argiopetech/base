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

#include "Star.hpp"
#include "densities.hpp"
#include "MsFilterSet.hpp"
#include "marg.hpp"
#include "Model.hpp"
#include "mpiMcmc.hpp"
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
vector<StellarSystem> readCmdData (struct ifmrMcmcControl &ctrl, const Model &evoModels, vector<int> &filters, std::array<double, FILTS> &filterPriorMin, std::array<double, FILTS> &filterPriorMax, const Settings &settings)
{
    string line, pch;

    vector<StellarSystem> systems;

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
                const_cast<Model&>(evoModels).numFilts++;
                break;
            }
        }
    }

    if (filters.empty())
    {
        cerr << "Exiting due to empty filter set. Did you specify the correct filterSet?" << endl;
        exit(1);
    }

    // This loop reads in photometry data
    // It also reads a best guess for the mass
    systems.clear();

    while (getline(ctrl.rData, line))
    {
        systems.emplace_back(line, filters.size());

        for (decltype(filters.size()) i = 0; i < filters.size(); ++i)
        {
            if (systems.back().obsPhot.at(i) < filterPriorMin.at(i))
            {
                filterPriorMin.at(i) = systems.back().obsPhot.at(i);
            }

            if (systems.back().obsPhot.at(i) > filterPriorMax.at(i))
            {
                filterPriorMax.at(i) = systems.back().obsPhot.at(i);
            }
        }

        if (!(systems.back().primary.status == 3 || (systems.back().obsPhot.at(settings.cluster.index) >= settings.cluster.minMag && systems.back().obsPhot.at(settings.cluster.index) <= settings.cluster.maxMag)))
        {
            systems.pop_back();
        }
    }

    return systems;
} /* readCmdData */


void printHeader (ofstream &file, array<double, NPARAMS> const &priors)
{
    const array<string, NPARAMS> paramNames = { "     logAge",
                                                "          Y",
                                                "        FeH",
                                                "    modulus",
                                                " absorption",
                                                "carbonicity",
                                                "  IFMRconst",
                                                "    IFMRlin",
                                                "   IFMRquad"};

    for (int p = 0; p < NPARAMS; p++)
    {
        if (priors.at(p) > EPSILON || p == MOD || p == FEH || p == ABS)
        {
            file << paramNames.at(p) << ' ';
        }
    }
    file << "logPost" << endl;
}

// Create Cholesky Decomp
void make_cholesky_decomp(struct ifmrMcmcControl &ctrl, MVatrix<double, NPARAMS> &params)
{
    double cov;
    int nParamsUsed = 0;

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
