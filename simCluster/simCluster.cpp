#include <array>
#include <vector>
#include <iostream>
#include <random>

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <boost/format.hpp>

#include "Cluster.hpp"
#include "Model.hpp"
#include "Settings.hpp"
#include "Star.hpp"
#include "FilterSet.hpp"

#include "WhiteDwarf.hpp"

using std::array;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

// Used by other methods in simCluster
static int nMSRG = 0, nWD = 0, nNSBH = 0;       //, nDa=0, nDb=0;
static double wdMassTotal = 0.0, MSRGMassTotal = 0.0;

double filterPriorMin[FILTS];
double filterPriorMax[FILTS];

// Used by a bunch of different functions.
vector<int> filters;

// For random # generator (mt19937ar.c)
unsigned long seed = 0;

int main (int argc, char *argv[])
{
    int i, nStars, cmpnt;//, nBrownDwarfs;
    double fractionBinary, tempU, massTotal, fractionDB, tempMod, minV, maxV, minMass = 0.15;
    char w_file[100];
    FILE *w_ptr;
    Cluster theCluster;
    StellarSystem theStar;

    double drawFromIMF (std::mt19937&);
    void updateCount (const StellarSystem&, int cmpnt);

    Settings settings;

    settings.fromCLI (argc, argv);
    if (!settings.files.config.empty())
    {
        settings.fromYaml (settings.files.config);
    }
    else
    {
        settings.fromYaml ("base9.yaml");
    }

    settings.fromCLI (argc, argv);

    const Model evoModels = makeModel(settings);

    for (int filt = 0; filt < 8; filt++)
        filters.push_back(filt); // Calculate all of U-K

//    theCluster.nStars = settings.simCluster.nStars;
    theCluster.setM_wd_up(settings.whiteDwarf.M_wd_up);
    fractionBinary = settings.simCluster.percentBinary;
    fractionDB = settings.simCluster.percentDB;
    theCluster.mod = settings.cluster.distMod;
    theCluster.abs = settings.cluster.Av;
    theCluster.age = settings.cluster.logClusAge;
    theCluster.feh = settings.cluster.Fe_H;
    theCluster.yyy = settings.cluster.Y;
    theCluster.carbonicity = settings.cluster.carbonicity;
//    nFieldStars = settings.simCluster.nFieldStars;
//    nBrownDwarfs = settings.simCluster.nBrownDwarfs;

    fractionBinary /= 100.;     // input as percentages, use as fractions
    fractionDB /= 100.;

    seed = settings.seed;

    // verbose = settings.verbose;
    // if (verbose < 0 || verbose > 2)
    //     verbose = 1;            // give standard feedback if incorrectly specified

    // !!! FIX ME !!!
    cerr << "This may be broken. If we need field stars and are using the YALE models, we also have to load the DSED models." << endl;
//    loadModels (theCluster, evoModels, settings);

    if (settings.mainSequence.msRgbModel == MsModel::YALE)
        minMass = 0.4;
    if (settings.mainSequence.msRgbModel == MsModel::DSED)
        minMass = 0.25;

    strcpy (w_file, settings.files.output.c_str());
    strcat (w_file, ".sim.out");
    if ((w_ptr = fopen (w_file, "w")) == NULL)
    {
        cerr << "\nFile " << w_file << " was not available for writing - exiting" << endl;
        exit (1);
    }

    nStars = 0;
    massTotal = 0.0;

    std::mt19937 gen(seed * uint32_t(2654435761));  // Applies Knuth's multiplicative hash for obfuscation (TAOCP Vol. 3)

    //Output headers
    fprintf (w_ptr, "id  mass1 ");
    for (auto f : filters)

        fprintf (w_ptr, "%s1 ", evoModels.filterSet->getFilterName(f).c_str());
    fprintf (w_ptr, "stage1 wdM1 wdType1 wdLogTeff1 ltau1 mass2 ");
    for (auto f : filters)
        fprintf (w_ptr, "%s2 ", evoModels.filterSet->getFilterName(f).c_str());
    fprintf (w_ptr, "stage2 wdM2 wdType2 wdLogTeff2 ltau2 ");
    for (auto f : filters)
        fprintf (w_ptr, "%s ", evoModels.filterSet->getFilterName(f).c_str());
    fprintf (w_ptr, "\n");

    minV = 1000.0;
    maxV = -1000.0;

    ////////////////////////////////
    ///// Create cluster stars /////
    ////////////////////////////////
    // derive masses, mags, and summary stats

    theCluster.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, theCluster.feh, theCluster.yyy, theCluster.age);    // determine AGBt ZAMS mass, to find evol state

    std::normal_distribution<double> normDist(1, 0.5);

    for (i = 0; i < settings.simCluster.nStars; i++)
    {                           // for all systems in the cluster
        do
        {
            theStar.primary.mass = drawFromIMF (gen); // create single stars in arbitrary mass order
        } while (theStar.primary.mass < minMass);
        nStars++;                       // keep track of the number of stars created
        massTotal += theStar.primary.mass;
        theStar.setMassRatio(0.0);

        fprintf (w_ptr, "%4d %7.3f ", i + 1, theStar.primary.mass); // output primary star data

        auto photometry = theStar.primary.getMags(theCluster, evoModels, filters);

        for (auto f : filters)
        {
                fprintf (w_ptr, "%6.3f ", photometry[f] < 99. ? photometry[f] : 99.999);
        }
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", theStar.primary.getStatus(theCluster), (theStar.primary.getStatus(theCluster) == WD ? theStar.primary.wdMassNow(theCluster, evoModels) : 0.0), 0, theStar.primary.wdLogTeff(theCluster, evoModels), theStar.primary.getLtau(theCluster, evoModels));

        tempU = theStar.primary.mass;              // save so we can output the total photometry below
        if (theStar.primary.getStatus(theCluster) != NSBH && theStar.primary.getStatus(theCluster) != WD)
        {
            if (std::generate_canonical<double, 53>(gen) < fractionBinary)
            {                           // create binaries among appropriate fraction
                do
                {
                    theStar.setMassRatio(normDist(gen));
                } while (theStar.getMassRatio() < 0 || theStar.getMassRatio() > 1);
                theStar.primary.mass = tempU * theStar.getMassRatio();
                nStars++;
                massTotal += theStar.primary.mass;
                theStar.setMassRatio(0.0);
            }
            else
                theStar.primary.mass = 0.0;
        }
        else
            theStar.primary.mass = 0.0;


        // Evolve secondary star by itself
        photometry = theStar.primary.getMags(theCluster, evoModels, filters);

        fprintf (w_ptr, "%7.3f ", theStar.primary.mass);    // output secondary star data
        for (auto f : filters)
        {
                fprintf (w_ptr, "%6.3f ", photometry[f] < 99. ? photometry[f] : 99.999);
        }
        // !! FIX ME !!
        // Should this really be the primary if we're outputting secondary star data?
        // There's no reason this can't use the secondary properly, now that a
        //   StellarSystem and a Star are seperate concepts.
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", theStar.primary.getStatus(theCluster), (theStar.primary.getStatus(theCluster) == WD ? theStar.primary.wdMassNow(theCluster, evoModels) : 0.0), 0, theStar.primary.wdLogTeff(theCluster, evoModels), theStar.primary.getLtau(theCluster, evoModels));

        theStar.secondary.mass = theStar.primary.mass;
        theStar.primary.mass = tempU;


        photometry = theStar.deriveCombinedMags(theCluster, evoModels, filters);

        switch (theStar.primary.getStatus(theCluster))
        {
            case MSRG:
                nMSRG++;
                MSRGMassTotal += theStar.primary.wdMassNow(theCluster, evoModels);
                break;
            case WD:
                nWD++;
                wdMassTotal += theStar.primary.wdMassNow(theCluster, evoModels);
                break;
            case NSBH:
                nNSBH++;
                break;
        }

        switch (theStar.secondary.getStatus(theCluster))
        {
            case MSRG:
                nMSRG++;
                MSRGMassTotal += theStar.secondary.wdMassNow(theCluster, evoModels);
                break;
            case WD:
                nWD++;
                wdMassTotal += theStar.secondary.wdMassNow(theCluster, evoModels);
                break;
            case NSBH:
                nNSBH++;
                break;
        }

        for (auto f : filters)
            fprintf (w_ptr, "%6.3f ", photometry[f] < 99. ? photometry[f] : 99.999);  // output photometry for whole system
        fprintf (w_ptr, "\n");

        // Update min and max
        if (photometry[2] < 97 && theStar.primary.getStatus(theCluster) != WD && theStar.secondary.getStatus(theCluster) != WD)
        {
            if (photometry[2] < minV)
                minV = photometry[2];
            if (photometry[2] > maxV)
                maxV = photometry[2];
        }
    }

    ///////////////////////////////
    ///// Create brown dwarfs /////
    ///////////////////////////////
/*
    for (i = 0; i < nBrownDwarfs; i++)
    {                           // for all systems in the cluster
        theStar.primary.mass = 0.0995 * std::generate_canonical<double, 53>(gen) + 0.0005; // create single stars in arbitrary mass order
        nStars++;                       // keep track of the number of stars created
        massTotal += theStar.primary.mass;
        theStar.massRatio = 0.0;
        theStar.primary.status = BD;

        evolve (theCluster, evoModels, globalMags, filters, theStar, ltau);      // given inputs, derive mags for first component

        fprintf (w_ptr, "%4d %7.4f ", i + 10001, theStar.primary.mass());     // output primary star data
        for (auto f : filters)
        {
                fprintf (w_ptr, "%6.3f ", theStar.photometry[f] < 99. ? theStar.photometry[f] : 99.999);
        }
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", theStar.primary.status, 0.0, 0, theStar.wdLogTeff[0], ltau[0]);

        fprintf (w_ptr, "%7.3f ", 0.00);        // output secondary star data
        for (auto f [[gnu::unused]] : filters)
        {
                fprintf (w_ptr, "%6.3f ", 99.999);
        }
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", DNE, 0.0, 0, 0.0, 0.0);

        for (cmpnt = 0; cmpnt < 2; cmpnt++)
            updateCount (&theStar, cmpnt);

        for (auto f : filters)
            fprintf (w_ptr, "%6.3f ", theStar.photometry[f] < 99. ? theStar.photometry[f] : 99.999);  // output photometry for whole system
        fprintf (w_ptr, "\n");
    }
*/
    cout << "\n Properties for cluster:" << endl;
    cout << boost::format(" logClusAge     = %6.3f") % theCluster.age << endl;
    cout << boost::format(" [Fe/H]         = %5.2f") % theCluster.feh << endl;
    cout << boost::format(" Y              = %5.2f") % theCluster.yyy << endl;
    cout << boost::format(" modulus        = %5.2f") % theCluster.mod << endl;
    cout << boost::format(" Av             = %5.2f") % theCluster.abs << endl;
    cout << boost::format(" WDMassUp       = %4.1f") % theCluster.getM_wd_up() << endl;
    cout << boost::format(" fractionBinary = %5.2f") % fractionBinary << endl;
    
    cout << "Totals:" << endl;
    cout << " nSystems       = " << settings.simCluster.nStars << endl;
    cout << " nStars         = " << nStars << endl;
    cout << " nMSRG          = " << nMSRG << endl;
    cout << " nWD            = " << nWD << endl;
    cout << " nNSBH          = " << nNSBH << endl;
    cout.precision(2);
    cout << " massTotal      = " << massTotal << endl;
    cout.precision(2);
    cout << " MSRGMassTotal  = " << MSRGMassTotal << endl;
    cout.precision(2);
    cout << " wdMassTotal    = " << wdMassTotal << endl;
    cout << " nFieldStars    = " << settings.simCluster.nFieldStars * 2 << endl;;


    //////////////////////////////
    ///// Create field stars /////
    //////////////////////////////

    double minFeH, maxFeH, minAge, maxAge;

//    theCluster.nStars = nFieldStars;
    tempMod = theCluster.mod;

    {
        // !!! FIX ME !!!
        cerr << "Field stars are broken..." << endl;
        exit(1);

        if (settings.mainSequence.msRgbModel == MsModel::YALE)
        {
            ; //            evoModels.mainSequenceEvol = DSED;
        }
    }

    if ((settings.mainSequence.msRgbModel == MsModel::DSED) || (settings.mainSequence.msRgbModel == MsModel::YALE))
    {
        minFeH = -2.5;
        maxFeH = 0.56;
        minAge = 8.4;
        maxAge = 10.17;
    }
    else
    {
        minFeH = -1.5;
        maxFeH = 0.2;
        minAge = 8.0;
        maxAge = 9.7;
    }

    std::vector<double> photometry;

    i = 0;
    while (i < settings.simCluster.nFieldStars)
    {
        do
        {
            do
            {
                theStar.primary.mass = drawFromIMF (gen);
            } while (theStar.primary.mass < minMass);
            theStar.setMassRatio(std::generate_canonical<double, 53>(gen));

            // Draw a new age and metallicity
            theCluster.age = minAge + (maxAge - minAge) * std::generate_canonical<double, 53>(gen);
            theCluster.feh = minFeH + (maxFeH - minFeH) * std::generate_canonical<double, 53>(gen);

            // Determine a new distance, weighted so
            // there are more stars behind than in front
            theCluster.mod = tempMod - 12.0 + log10 (exp10 ((pow (pow (26.0, 3.0) * std::generate_canonical<double, 53>(gen), 1.0 / 3.0))));

            theCluster.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, theCluster.feh, theCluster.yyy, theCluster.age);    // determine AGBt ZAMS mass, to find evol state

            photometry = theStar.deriveCombinedMags(theCluster, evoModels, filters);

        } while (photometry[2] < minV || photometry[2] > maxV || photometry[1] - photometry[2] < -0.5 || photometry[1] - photometry[2] > 1.7);

        fprintf (w_ptr, "%4d %7.3f ", i + 20001, theStar.primary.mass);
        for (auto f [[gnu::unused]] : filters)
            fprintf (w_ptr, "%6.3f ", 99.999);
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", theStar.primary.getStatus(theCluster), (theStar.primary.getStatus(theCluster) == WD ? theStar.primary.wdMassNow(theCluster, evoModels) : 0.0), 0, theStar.primary.wdLogTeff(theCluster, evoModels), theStar.primary.getLtau(theCluster, evoModels));
        fprintf (w_ptr, "%7.3f ", theStar.secondary.mass);
        for (auto f [[gnu::unused]] : filters)
            fprintf (w_ptr, "%6.3f ", 99.999);
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", theStar.secondary.getStatus(theCluster), 0.0, 0, theStar.secondary.wdLogTeff(theCluster, evoModels), theStar.secondary.getLtau(theCluster, evoModels));
        for (auto f : filters)
            fprintf (w_ptr, "%9.6f ", photometry[f]);
        fprintf (w_ptr, "\n");

        i++;
    }

    fclose (w_ptr);

    return (0);
}

void updateCount (const StellarSystem &system, int cmpnt)
{
    Star star;

    if (cmpnt == 0)
        star = system.primary;
    else
        star = system.secondary;


}
