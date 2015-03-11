#include <Rcpp.h>
#include <iostream>
#include <string>
#include <vector>

#include "Cluster.hpp"
#include "constants.hpp"
#include "marg.hpp"
#include "Matrix.hpp"
#include "Model.hpp"
#include "ifmr.hpp"
#include "Star.hpp"

namespace evil {
    Settings s;
    bool initialized = false;

    class globals
    {
      private:
        globals()
            : evoModels(makeModel(s))
        {
            evil::initialized = true;

            for (auto &a : clust.priorVar)
                a = std::numeric_limits<double>::max();
            for (auto &a : clust.priorMean)
                a = std::numeric_limits<double>::max();
            for (auto &a : clust.mean)
                a = std::numeric_limits<double>::max();

            std::vector<std::string> msFilters = evoModels.mainSequenceEvol->getAvailableFilters();
            std::vector<std::string> wdFilters = evoModels.WDAtmosphere->getAvailableFilters();

            {
                for ( auto m : msFilters )
                {
                    for ( auto w : wdFilters )
                    {
                        if (m == w)
                        {
                            filters.push_back(m);
                            break;
                        }
                    }
                }
            }

            for (size_t f = 0; f < filters.size(); ++f)
            {
                try
                {
                    Filters::absCoeffs.at(filters.at(f));
                }
                catch(std::out_of_range &e)
                {
                    filters.erase(filters.begin() + f);

                    if (f > 0)
                        f -= 1;
                }
            }

            evoModels.restrictFilters(filters);

//            clust.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(clust.feh, clust.yyy, clust.age);    // determine AGBt ZAMS mass, to find evol state
        }
        globals(const globals&) = delete;
        globals(const globals&&) = delete;
        globals& operator=(const globals&) = delete;

      public:
        Cluster clust;
        Model evoModels;
        std::vector<std::string> filters;

        std::unique_ptr<Isochrone> isochrone;

        static globals& getInstance()
        {
            static globals instance;

            return instance;
        }

        void reloadModels()
        {
            evoModels = makeModel(s);
        }
    };
}

const double LOG_G_PLUS_LOG_M_SUN = 26.12302173752;

bool isInitialized()
{
    if (evil::initialized)
    {
        return true;
    }
    else
    {
        throw Rcpp::exception("You need to initialize first...");
        return false; // This will never get called...
    }
}


// [[Rcpp::export]]
std::vector<std::string> listFilters()
{
    return evil::globals::getInstance().filters;
}


// [[Rcpp::export]]
double getAGBt_zmass()
{
    return evil::globals::getInstance().isochrone->agbTipMass();
}


// [[Rcpp::export]]
void initBase(std::string modelDir, int msModel, int wdModel, int ifmr)
{
    evil::s.files.models = modelDir;

    evil::s.mainSequence.msRgbModel = static_cast<MsModel>(msModel);

    evil::s.whiteDwarf.ifmr = ifmr;
    evil::s.whiteDwarf.wdModel = static_cast<WdModel>(wdModel);
    evil::s.whiteDwarf.M_wd_up = 8.0;

    evil::s.cluster.priorMeans.Fe_H = 0.07;
    evil::s.cluster.priorSigma.Fe_H = 0.05;

    evil::s.cluster.priorMeans.distMod = 0.0;
    evil::s.cluster.priorSigma.distMod = 0.05;

    evil::s.cluster.priorMeans.Av = 0.009;
    evil::s.cluster.priorSigma.Av = 0.006;

    evil::s.cluster.priorMeans.Y = 0.29;
    evil::s.cluster.priorSigma.Y = 0.0;

    evil::s.cluster.priorMeans.carbonicity = 0.5;
    evil::s.cluster.priorSigma.carbonicity = 0.038;

    evil::s.cluster.priorMeans.logAge = 8.796;

    evil::s.cluster.minMag = 0.0;
    evil::s.cluster.maxMag = 30.0;
    evil::s.cluster.index = 2;

    evil::globals::getInstance();
}


// [[Rcpp::export]]
void setClusterParameters(double age, double feh, double distMod, double av, double y, double carbonicity)
{
    if (isInitialized())
    {
        const auto evoModels = evil::globals::getInstance().evoModels;
        auto &clust = evil::globals::getInstance().clust;
        auto &isochrone = evil::globals::getInstance().isochrone;

        clust.age = age;
        clust.feh = feh;
        clust.mod = distMod;
        clust.abs = av;
        clust.yyy = y;
        clust.carbonicity = carbonicity;

        isochrone.reset(evoModels.mainSequenceEvol->deriveIsochrone(clust.feh, clust.yyy, clust.age));
    }
}

// [[Rcpp::export]]
void changeModels(int msModel, int wdModel, int ifmr)
{
    if (isInitialized())
    {
        evil::s.mainSequence.msRgbModel = static_cast<MsModel>(msModel);

        evil::s.whiteDwarf.ifmr = ifmr;
        evil::s.whiteDwarf.wdModel = static_cast<WdModel>(wdModel);
        evil::s.whiteDwarf.M_wd_up = 8.0;

        const auto evoModels = evil::globals::getInstance().evoModels;
        auto &clust = evil::globals::getInstance().clust;
        auto &isochrone = evil::globals::getInstance().isochrone;

        evil::globals::getInstance().reloadModels();

        isochrone.reset(evoModels.mainSequenceEvol->deriveIsochrone(clust.feh, clust.yyy, clust.age));
    }
}


// [[Rcpp::export]]
void setIFMRParameters(double intercept, double slope, double quadCoef)
{
    if (isInitialized())
    {
        evil::globals::getInstance().clust.ifmrIntercept = intercept;
        evil::globals::getInstance().clust.ifmrSlope = slope;
        evil::globals::getInstance().clust.ifmrQuadCoef = quadCoef;
    }
}


// [[Rcpp::export]]
std::vector<double> evolve (double mass1, double mass2)
{
    if (isInitialized())
    {
        double flux;
        StellarSystem system;
        system.primary.mass = mass1;
        system.secondary.mass = mass2;

        const auto evoModels = evil::globals::getInstance().evoModels;
        const auto clust = evil::globals::getInstance().clust;
        auto &isochrone = evil::globals::getInstance().isochrone;

        // AGBt_zmass never set because age and/or metallicity out of range of models.
        if (isochrone->agbTipMass() < EPS)
        {
            throw Rcpp::exception("Bounds error in evolve");
        }

        return system.deriveCombinedMags(clust, evoModels, *isochrone);
    }
    else
    {
        throw Rcpp::exception("This should never happen...");
    }
}
