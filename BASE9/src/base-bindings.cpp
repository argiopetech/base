#include <Rcpp.h>
#include <array>
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

            clust.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, clust.feh, clust.yyy, clust.age);    // determine AGBt ZAMS mass, to find evol state
        }
        globals(const globals&) = delete;
        globals(const globals&&) = delete;
        globals& operator=(const globals&) = delete;

      public:
        Cluster clust;
        Model evoModels;
        std::vector<int> filters = {0, 1, 2, 3, 4, 5, 6, 7};

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
void initBase(std::string modelDir, int msFilter, int msModel, int wdModel, int ifmr)
{
    evil::s.files.models = modelDir;

    evil::s.mainSequence.filterSet = static_cast<MsFilter>(msFilter);
    evil::s.mainSequence.msRgbModel = static_cast<MsModel>(msModel);

    evil::s.whiteDwarf.ifmr = ifmr;
    evil::s.whiteDwarf.wdModel = static_cast<WdModel>(wdModel);
    evil::s.whiteDwarf.M_wd_up = 8.0;

    evil::s.cluster.Fe_H = 0.07;
    evil::s.cluster.sigma.Fe_H = 0.05;

    evil::s.cluster.distMod = 0.0;
    evil::s.cluster.sigma.distMod = 0.05;

    evil::s.cluster.Av = 0.009;
    evil::s.cluster.sigma.Av = 0.006;

    evil::s.cluster.Y = 0.29;
    evil::s.cluster.sigma.Y = 0.0;

    evil::s.cluster.carbonicity = 0.5;
    evil::s.cluster.sigma.carbonicity = 0.038;

    evil::s.cluster.logClusAge = 8.796;

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
        const auto filters = evil::globals::getInstance().filters;
        const auto evoModels = evil::globals::getInstance().evoModels;
        auto &clust = evil::globals::getInstance().clust;

        clust.age = age;
        clust.feh = feh;
        clust.mod = distMod;
        clust.abs = av;
        clust.yyy = y;
        clust.carbonicity = carbonicity;

        clust.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, clust.feh, clust.yyy, clust.age);    // determine AGBt ZAMS mass, to find evol state
    }
}

// [[Rcpp::export]]
void changeModels(int msFilter, int msModel, int wdModel, int ifmr)
{
    if (isInitialized())
    {
        evil::s.mainSequence.filterSet = static_cast<MsFilter>(msFilter);
        evil::s.mainSequence.msRgbModel = static_cast<MsModel>(msModel);

        evil::s.whiteDwarf.ifmr = ifmr;
        evil::s.whiteDwarf.wdModel = static_cast<WdModel>(wdModel);
        evil::s.whiteDwarf.M_wd_up = 8.0;

        const auto filters = evil::globals::getInstance().filters;
        const auto evoModels = evil::globals::getInstance().evoModels;
        auto &clust = evil::globals::getInstance().clust;

        evil::globals::getInstance().reloadModels();

        clust.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, clust.feh, clust.yyy, clust.age);    // determine AGBt ZAMS mass, to find evol state
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
std::array<double, 8> evolve (double mass1, double mass2)
{
    if (isInitialized())
    {
        std::array<double, 8> returnMags;
        Matrix<double, 2, 14> mags;
        std::array<double, 14> combinedMags;

        double flux;
        StellarSystem system;
        system.primary.mass = mass1;
        system.secondary.mass = mass2;

        for (auto &mag : mags)
        {
            mag.fill(0.0);
        }

        combinedMags.fill(0.0);

        const auto filters = evil::globals::getInstance().filters;
        const auto evoModels = evil::globals::getInstance().evoModels;
        const auto clust = evil::globals::getInstance().clust;

        // AGBt_zmass never set because age and/or metallicity out of range of models.
        if (clust.AGBt_zmass < EPS)
        {
            throw Rcpp::exception("Bounds error in evolve");
        }

        combinedMags = system.deriveCombinedMags (clust, evoModels, filters);

        for (auto f : filters)
        {
            returnMags[f] = combinedMags[f];
        }

        return returnMags;
    }
    else
    {
        throw Rcpp::exception("This should never happen...");
    }
}
