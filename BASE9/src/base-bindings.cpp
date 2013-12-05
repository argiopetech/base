#include <Rcpp.h>
#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "Cluster.hpp"
#include "constants.hpp"
#include "gBergMag.hpp"
#include "Matrix.hpp"
#include "Model.hpp"
#include "ifmr.hpp"

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

            for (auto &a : pCluster.betaAgeMod)
                a = std::numeric_limits<double>::max();
            for (auto &a : pCluster.priorVar)
                a = std::numeric_limits<double>::max();
            for (auto &a : pCluster.priorMean)
                a = std::numeric_limits<double>::max();
            for (auto &a : pCluster.mean)
                a = std::numeric_limits<double>::max();

            pCluster.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, pCluster.feh, pCluster.yyy, pCluster.age);    // determine AGBt ZAMS mass, to find evol state
        }
        globals(const globals&) = delete;
        globals(const globals&&) = delete;
        globals& operator=(const globals&) = delete;

      public:
        Cluster pCluster;
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
    const auto pCluster = evil::globals::getInstance().pCluster;

    std::pair<double, double> teffRadiusPair;

    double thisWDMass = 0.0, thisPrecLogAge = 0.0, thisLogTeff, thisWDLogRadius = 0.0, thisWDLogG = 0.0;

    thisPrecLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(pCluster.feh, mass);

    thisWDMass = intlFinalMassReln (pCluster, evoModels, mass);

    //get temperature from WD cooling models (returns 0.0 if there is an error(or does it??))
    teffRadiusPair = evoModels.WDcooling->wdMassToTeffAndRadius (pCluster.age, pCluster.carbonicity, thisPrecLogAge, thisWDMass);

    thisLogTeff = teffRadiusPair.first;
    thisWDLogRadius = teffRadiusPair.second;

    //*******this now gets trapped for in wdMassToTeffAndRadius so it should be unnecessary here (???)
    if (thisPrecLogAge >= pCluster.age)
    {                           // mcmc.c can cause this by adjusting masses and ages         
        for (auto f : filters)
            globalMags[f] = -4.; // place at tip of RGB                                       
    }
    else
    {
        //Calculate logg                                                                      
        thisWDLogG = LOG_G_PLUS_LOG_M_SUN + log10 (thisWDMass) - 2 * thisWDLogRadius;
        bergeronTeffToMags (filters, globalMags, thisLogTeff, thisWDLogG, DA); // All WD stars are currently set to DA
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
        auto &pCluster = evil::globals::getInstance().pCluster;

        evil::globals::getInstance().pCluster.age = age;
        evil::globals::getInstance().pCluster.feh = feh;
        evil::globals::getInstance().pCluster.mod = distMod;
        evil::globals::getInstance().pCluster.abs = av;
        evil::globals::getInstance().pCluster.yyy = y;
        evil::globals::getInstance().pCluster.carbonicity = carbonicity;

        pCluster.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, pCluster.feh, pCluster.yyy, pCluster.age);    // determine AGBt ZAMS mass, to find evol state
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
        auto &pCluster = evil::globals::getInstance().pCluster;

        evil::globals::getInstance().reloadModels();

        pCluster.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, pCluster.feh, pCluster.yyy, pCluster.age);    // determine AGBt ZAMS mass, to find evol state
    }
}


// [[Rcpp::export]]
void setIFMRParameters(double intercept, double slope, double quadCoef)
{
    if (isInitialized())
    {
        evil::globals::getInstance().pCluster.ifmrIntercept = intercept;
        evil::globals::getInstance().pCluster.ifmrSlope = slope;
        evil::globals::getInstance().pCluster.ifmrQuadCoef = quadCoef;
    }
}


// [[Rcpp::export]]
std::array<double, 8> evolve (double mass)
{
    if (isInitialized())
    {
        std::array<double, 14> globalMags;
        std::array<double, 8> mag;

        const auto filters = evil::globals::getInstance().filters;
        const auto evoModels = evil::globals::getInstance().evoModels;
        const auto pCluster = evil::globals::getInstance().pCluster;

        // AGBt_zmass never set because age and/or metallicity out of range of models.
        if (pCluster.AGBt_zmass < EPS)
        {
            throw Rcpp::exception("Bounds error in evolve");
        }

        if (mass <= 0.0001)
        {                           // for non-existent secondary stars
            for (auto f : filters)
                mag[f] = 99.999;
        }
        else if (mass <= pCluster.AGBt_zmass)
        {                           // for main seq or giant star
            evoModels.mainSequenceEvol->msRgbEvol(filters, globalMags, mass);
            for (auto f : filters)
                mag[f] = globalMags[f];
        }
        else if (mass <= pCluster.M_wd_up)
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

        return mag;
    }
    else
    {
        throw Rcpp::exception("This should never happen...");
    }
}
