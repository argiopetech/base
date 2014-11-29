#include <array>
#include <chrono>
#include <vector>
#include <iostream>
#include <iomanip>
#include <memory>
#include <string>
#include <random>

#include <cstring>

#include "Cluster.hpp"
#include "Model.hpp"
#include "Settings.hpp"
#include "Star.hpp"
#include "Filters.hpp"

#include "WhiteDwarf.hpp"

using std::array;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::fixed;
using std::setprecision;
using std::setw;
using std::string;
using std::unique_ptr;

// Used by other methods in simCluster
static int nMSRG = 0, nWD = 0, nNSBH = 0;       //, nDa=0, nDb=0;
static double wdMassTotal = 0.0, MSRGMassTotal = 0.0;

// For random # generator (mt19937ar.c)
unsigned long seed = 0;

int main (int argc, char *argv[])
{
    int i, nStars;//, nBrownDwarfs;
    double fractionBinary, massTotal, fractionDB, tempMod, minV, maxV, minMass = 0.15;
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

    if (settings.seed == std::numeric_limits<uint32_t>::max())
    {
        srand(std::time(0));
        settings.seed = rand();

        cout << "Seed: " << settings.seed << endl;
    }

    Model evoModels = makeModel(settings);

    vector<string> filters;

    {
        vector<string> msFilters = evoModels.mainSequenceEvol->getAvailableFilters();
        vector<string> wdFilters = evoModels.WDAtmosphere->getAvailableFilters();

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

    theCluster.setM_wd_up(settings.whiteDwarf.M_wd_up);
    fractionBinary = settings.simCluster.percentBinary;
    fractionDB = settings.simCluster.percentDB;
    theCluster.mod = settings.cluster.distMod;
    theCluster.abs = settings.cluster.Av;
    theCluster.age = settings.cluster.logAge;
    theCluster.feh = settings.cluster.Fe_H;
    theCluster.yyy = settings.cluster.Y;
    theCluster.carbonicity = settings.cluster.carbonicity;

    fractionBinary /= 100.;     // input as percentages, use as fractions
    fractionDB /= 100.;

    seed = settings.seed;

    if (settings.mainSequence.msRgbModel == MsModel::YALE)
        minMass = 0.4;
    if (settings.mainSequence.msRgbModel == MsModel::OLD_DSED)
        minMass = 0.25;

    strcpy (w_file, (settings.files.output + ".sim.out").c_str());
    if ((w_ptr = fopen (w_file, "w")) == NULL)
    {
        cerr << "\nFile " << w_file << " was not available for writing - exiting" << endl;
        exit (1);
    }

    nStars = 0;
    massTotal = 0.0;

    std::mt19937 gen(seed * uint32_t(2654435761));  // Applies Knuth's multiplicative hash for obfuscation (TAOCP Vol. 3)

    //Output headers
    // Star     U       B       V      sigU    sigB    sigV    mass1 massRatio stage Cmprior useDBI
    fprintf (w_ptr, "%4s ", "id");

    for (auto f : filters)
        fprintf (w_ptr, "%6s ", f.c_str());

    for (auto f : filters)
        fprintf (w_ptr, "%6s ", ("sig" + f).c_str());

    fprintf(w_ptr, "%7s %9s %5s %7s %6s\n", "mass1", "massRatio", "stage", "Cmprior", "useDBI");

    minV = 1000.0;
    maxV = -1000.0;

    ////////////////////////////////
    ///// Create cluster stars /////
    ////////////////////////////////
    // derive masses, mags, and summary stats

    unique_ptr<Isochrone> isochrone(evoModels.mainSequenceEvol->deriveIsochrone(theCluster.feh, theCluster.yyy, theCluster.age));

    std::normal_distribution<double> normDist(1, 0.5);

    for (i = 0; i < settings.simCluster.nStars; i++)
    {                           // for all systems in the cluster
        do
        {
            theStar.primary.mass = drawFromIMF (gen); // create single stars in arbitrary mass order
        } while (theStar.primary.mass < minMass);
        nStars++;                       // keep track of the number of stars created
        massTotal += theStar.primary.mass;

        if (theStar.primary.getStatus(theCluster, *isochrone) != NSBH && theStar.primary.getStatus(theCluster, *isochrone) != WD)
        {
            if (std::generate_canonical<double, 53>(gen) < fractionBinary)
            {                           // create binaries among appropriate fraction
                double secondaryMass;
                do
                {
                    secondaryMass = normDist(gen);
                } while (secondaryMass < 0 || secondaryMass > 1);

                theStar.setMassRatio(secondaryMass);

                nStars++;
                massTotal += theStar.secondary.mass;
            }
            else
                theStar.secondary.mass = 0.0;
        }
        else
            theStar.secondary.mass = 0.0;

        fprintf (w_ptr, "%4d ", i + 1); // output primary star data

        // Evolve secondary star by itself
        auto photometry = theStar.deriveCombinedMags(theCluster, evoModels, *isochrone);

        for (size_t f = 0; f < filters.size(); ++f)
            fprintf (w_ptr, "%6.3f ", photometry[f] < 99. ? photometry[f] : 99.999);  // output photometry for whole system

        for (size_t f = 0; f < filters.size(); ++f)
            fprintf (w_ptr, "%6.3f ", 0.0);  // output (blank) sigmas

        fprintf (w_ptr, "%7.3f %9.3f %5d %7.3f %6d", theStar.primary.mass, theStar.getMassRatio(), theStar.primary.getStatus(theCluster, *isochrone), 0.999, 1);    // output secondary star data


        switch (theStar.primary.getStatus(theCluster, *isochrone))
        {
            case MSRG:
                nMSRG++;
                MSRGMassTotal += theStar.primary.wdMassNow(theCluster, evoModels, *isochrone);
                break;
            case WD:
                nWD++;
                wdMassTotal += theStar.primary.wdMassNow(theCluster, evoModels, *isochrone);
                break;
            case NSBH:
                nNSBH++;
                break;
        }

        switch (theStar.secondary.getStatus(theCluster, *isochrone))
        {
            case MSRG:
                nMSRG++;
                MSRGMassTotal += theStar.secondary.wdMassNow(theCluster, evoModels, *isochrone);
                break;
            case WD:
                nWD++;
                wdMassTotal += theStar.secondary.wdMassNow(theCluster, evoModels, *isochrone);
                break;
            case NSBH:
                nNSBH++;
                break;
        }


        fprintf (w_ptr, "\n");

        // Update min and max
        if (photometry[2] < 97 && theStar.primary.getStatus(theCluster, *isochrone) != WD && theStar.secondary.getStatus(theCluster, *isochrone) != WD)
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
    cout << "\nProperties for cluster:" << endl;
    cout << " logAge         = " << setw(6) << setprecision(3) << fixed << theCluster.age << endl;
    cout << " [Fe/H]         = " << setw(5) << setprecision(2) << fixed << theCluster.feh << endl;
    cout << " Y              = " << setw(5) << setprecision(2) << fixed << theCluster.yyy << endl;
    cout << " modulus        = " << setw(5) << setprecision(2) << fixed << theCluster.mod << endl;
    cout << " Av             = " << setw(5) << setprecision(2) << fixed << theCluster.abs << endl;
    cout << " WDMassUp       = " << setw(4) << setprecision(1) << fixed << theCluster.getM_wd_up() << endl;
    cout << " fractionBinary = " << setw(5) << setprecision(2) << fixed << fractionBinary << endl;
    
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
/*    cout << " nFieldStars    = " << settings.simCluster.nFieldStars * 2 << endl;;


    //////////////////////////////
    ///// Create field stars /////
    //////////////////////////////

    double minFeH, maxFeH, minAge, maxAge;

//    theCluster.nStars = nFieldStars;
    tempMod = theCluster.mod;

    if ((settings.mainSequence.msRgbModel == MsModel::OLD_DSED) || (settings.mainSequence.msRgbModel == MsModel::YALE))
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

            isochrone.reset(evoModels.mainSequenceEvol->deriveIsochrone(theCluster.feh, theCluster.yyy, theCluster.age));

            photometry = theStar.deriveCombinedMags(theCluster, evoModels, *isochrone);

        } while (photometry[2] < minV || photometry[2] > maxV || photometry[1] - photometry[2] < -0.5 || photometry[1] - photometry[2] > 1.7);

        fprintf (w_ptr, "%4d %7.3f ", i + 20001, theStar.primary.mass);
        for (auto f [[gnu::unused]] : filters)
            fprintf (w_ptr, "%6.3f ", 99.999);
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", theStar.primary.getStatus(theCluster, *isochrone), (theStar.primary.getStatus(theCluster, *isochrone) == WD ? theStar.primary.wdMassNow(theCluster, evoModels, *isochrone) : 0.0), 0, theStar.primary.wdLogTeff(theCluster, evoModels), theStar.primary.getLtau(theCluster, evoModels));
        fprintf (w_ptr, "%7.3f ", theStar.secondary.mass);
        for (auto f [[gnu::unused]] : filters)
            fprintf (w_ptr, "%6.3f ", 99.999);
        fprintf (w_ptr, "%d %5.3f %d %5.3f %5.3f ", theStar.secondary.getStatus(theCluster, *isochrone), 0.0, 0, theStar.secondary.wdLogTeff(theCluster, evoModels), theStar.secondary.getLtau(theCluster, evoModels));
        for (size_t f = 0; f < filters.size(); ++f)
            fprintf (w_ptr, "%9.6f ", photometry[f]);
        fprintf (w_ptr, "\n");

        i++;
    }
*/
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
