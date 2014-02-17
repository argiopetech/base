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

            for (auto &a : clust.betaAgeMod)
                a = std::numeric_limits<double>::max();
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

double wdEvol (std::array<double, 14> &globalMags, double mass)
{
    const auto filters = evil::globals::getInstance().filters;
    const auto evoModels = evil::globals::getInstance().evoModels;
    const auto clust = evil::globals::getInstance().clust;

    std::pair<double, double> teffRadiusPair;

    double thisWDMass = 0.0, thisPrecLogAge = 0.0, thisLogTeff, thisWDLogRadius = 0.0, thisWDLogG = 0.0;

    thisPrecLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(clust.feh, mass);

    thisWDMass = intlFinalMassReln (clust, evoModels, mass);

    //get temperature from WD cooling models (returns 0.0 if there is an error(or does it??))
    teffRadiusPair = evoModels.WDcooling->wdMassToTeffAndRadius (clust.age, clust.carbonicity, thisPrecLogAge, thisWDMass);

    thisLogTeff = teffRadiusPair.first;
    thisWDLogRadius = teffRadiusPair.second;

    //*******this now gets trapped for in wdMassToTeffAndRadius so it should be unnecessary here (???)
    if (thisPrecLogAge >= clust.age)
    {                           // mcmc.c can cause this by adjusting masses and ages         
        for (auto f : filters)
            globalMags[f] = -4.; // place at tip of RGB                                       
    }
    else
    {
        //Calculate mass
        globalMags = evoModels.WDAtmosphere->teffToMags (thisLogTeff, thisWDMass, WdAtmosphere::DA); // All WD stars are currently set to DA
    }

    return thisPrecLogAge;
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

        evil::globals::getInstance().clust.age = age;
        evil::globals::getInstance().clust.feh = feh;
        evil::globals::getInstance().clust.mod = distMod;
        evil::globals::getInstance().clust.abs = av;
        evil::globals::getInstance().clust.yyy = y;
        evil::globals::getInstance().clust.carbonicity = carbonicity;

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
        std::array<double, 8> combinedMags;
        Matrix<double, 3, 14> mags;

        double flux;
        Star pStar;

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

        auto mySetMags = [=](std::array<double, 14> &mag, double mass)
        {
            std::array<double, 14> globalMags;
            globalMags.fill(0.0);

            if (mass <= 0.0001)
            {                           // for non-existent secondary stars
                for (auto f : filters)
                    mag[f] = 99.999;
            }
            else if (mass <= clust.AGBt_zmass)
            {                           // for main seq or giant star
                evoModels.mainSequenceEvol->msRgbEvol(filters, globalMags, mass);
                for (auto f : filters)
                    mag[f] = globalMags[f];
            }
            else if (mass <= clust.M_wd_up)
            {                           // for white dwarf
                wdEvol (globalMags, mass);
                for (auto f : filters)
                    mag[f] = globalMags[f];
            }
            else if (mass <= 100.)
            {                           // for neutron star or black hole remnant
                for (auto f : filters)
                    mag[f] = 99.999;
            }
            else
            { // This should never happen...
                for (auto f : filters)
                    mag[f] = 99.999;
            }
        };

        mySetMags(mags[0], mass1);
        mySetMags(mags[1], mass2);

        deriveCombinedMags(mags, clust.abs, flux, clust, pStar, evoModels, filters);

        for (auto f : filters)
        {
            combinedMags[f] = mags[2][f];
        }

        return combinedMags;
    }
    else
    {
        throw Rcpp::exception("This should never happen...");
    }
}
