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
    stars.clear();

    while (getline(ctrl.rData, line))
    {
        stars.push_back(Star(line, filters.size()));

        for (decltype(filters.size()) i = 0; i < filters.size(); ++i)
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
        mc.clust.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, mc.clust.feh, mc.clust.yyy, mc.clust.age);    // determine AGBt ZAMS mass, to find evol state
        evolve (mc.clust, evoModels, globalMags, filters, star, ltau);

        if (star.status[0] == WD)
        {
            star.UStepSize = 0.05;      // use larger initial step size for white dwarfs
            star.massRatio = 0.0;
        }
    }
} /* initChain */


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
