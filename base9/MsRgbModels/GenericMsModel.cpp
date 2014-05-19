#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "Cluster.hpp"
#include "Isochrone.hpp"
#include "LinearTransform.hpp"
#include "Matrix.hpp"
#include "GenericMsModel.hpp"

using std::array;
using std::cerr;
using std::endl;
using std::ifstream;
using std::lower_bound;
using std::string;
using std::stringstream;
using std::vector;

const unsigned int maxIgnore = std::numeric_limits<char>::max();

void GenericMsModel::restrictToFilters(const vector<string>& filters)
{
    vector<int> indices;

    for (auto f : filters)
    {
        bool foundFilter = false;

        for (size_t i = 0; i < availableFilters.size(); ++i)
        {
            if ( f == availableFilters.at(i) )
            {
                indices.push_back(i);

                foundFilter = true;
                break;
            }
        }

        if ( ! foundFilter )
        {
            cout << "Couldn't find filter \"" << f << "\" in selected model set" << endl;
            exit(-1);
        }
    }

    for (size_t i = 0; i < indices.size(); ++i)
    {
        // Prepend the desired filter values to the beginning of the list...
        availableFilters.insert(availableFilters.begin() + i, availableFilters.at(indices.at(i) + i));
    }

    // ...then erase the rest.
    availableFilters.erase(availableFilters.begin() + indices.size(), availableFilters.end());
    // And try to reduce the capacity() to what is actually needed
    availableFilters.shrink_to_fit();

    // We should now have the same number of filters in the available filter set as in the input vector.
    assert(availableFilters.size() == filters.size());

    // Having this as a seperate loop allows a size assertion without
    // going over the data structure again. It also probably isn't bad
    // for cache locality (not that it matters)
    for (auto &fehs : fehCurves)
    {
        for (auto &ys : fehs.heliumCurves)
        {
            for (auto &is : ys.isochrones)
            {
                for (auto &eep : is.eeps)
                {
                    for (size_t i = 0; i < indices.size(); ++i)
                    {
                        // Prepend the desired filter values to the beginning of the list...
                        eep.mags.insert(eep.mags.begin() + i, eep.mags.at(indices.at(i) + i));
                    }

                    // ...then erase the rest.
                    eep.mags.erase(eep.mags.begin() + indices.size(), eep.mags.end());
                    // And try to reduce the capacity() to what is actually needed
                    eep.mags.shrink_to_fit();

                    // At this point, eep.mags, availableFilters, and filters should all be the same size.
                    assert(eep.mags.size() <= filters.size());
                }
            }
        }
    }
}


void GenericMsModel::loadModel(string path, FilterSetName)
{
    ifstream fin;

    bool haveFirstRecord = false;

    double logAge
        , previousFeH, currentFeH
        , previousY, currentY
        , previousAlphaFe, currentAlphaFe
        , previousMixLen, currentMixLen;

    string line;
    stringstream lin;

    vector<HeliumCurve> heliumCurves;
    vector<Isochrone> isochrones;
    vector<EvolutionaryPoint> eeps;

    string file = getFileName(path);

    // Open the file
    fin.open(file);

    // Ensure the file opened successfully
    if (!fin)
    {
        cerr << "\n file " << file << " was not found - exiting" << endl;
        exit (1);
    }

    // Now, for every line in the file...
    while (!fin.eof())
    {
        getline(fin, line);

        // Ignore comments (lines starting with '#')
        // The line should never be empty, but we segfault
        //   if it is (and we don't check), so it's better to check.
        if (!line.empty() && (line.at(0) != '#'))
        {
            stringstream in(line);

            if (line.at(0) == '%') // It's a command line
            {
                // Push the previous age
                if (haveFirstRecord)
                {
                    isochrones.emplace_back(logAge, eeps);
                    eeps.clear();
                }

                // A few choices...
                if (line.at(1) == 'f') // The once-per-file filter list
                {
                    string s;
                    in >> s; // Eat the command string "%f"

                    assert(s == "%f");

                    // While the stream is still good
                    while (in)
                    {
                        in >> s; // Input a filter
                        availableFilters.push_back(s); // And declare it "available"
                        in.peek(); // Try the next input to trigger end of stream
                        // This keeps us from getting a duplicate string at the end of the vector.
                    }
                }
                else if (line.at(1) == 's') // The section headers
                {
                    // %s [Fe/H]= 0.000000    [alpha/Fe]=0.000000    l/Hp=1.938000    Y=0.270000
                    in.ignore(maxIgnore, '='); // Kill everything up to the [Fe/H] value
                    in >> currentFeH;

                    in.ignore(maxIgnore, '='); // Kill everything up to the [alpha/Fe] value
                    in >> currentAlphaFe;

                    in.ignore(maxIgnore, '='); // Kill everything up to the l/Hp value
                    in >> currentMixLen;

                    in.ignore(maxIgnore, '='); // Kill everything up to the Y value
                    in >> currentY;

                    // Push the previous section
                    if (haveFirstRecord)
                    {
                        // We've hit a new FeH
                        if (currentFeH != previousFeH)
                        {
                            heliumCurves.emplace_back(previousY, isochrones);
                            fehCurves.emplace_back(previousFeH, heliumCurves);

                            // Clear everything
                            heliumCurves.clear();
                            isochrones.clear();
                            eeps.clear();
                            haveFirstRecord = false;
                        }
                        else if (currentY != previousY)
                        {
                            heliumCurves.emplace_back(previousY, isochrones);

                            // Clear everything except heliumCurves
                            isochrones.clear();
                            eeps.clear();
                        }
                        else // Something changed that we don't support
                        {
                            cerr << "Something changed that we don't support:\n\t" << line << endl;
                            exit(1);
                        }
                    }

                    // Save all the values for the next section header check
                    previousFeH = currentFeH;
                    previousAlphaFe = currentAlphaFe;
                    previousMixLen = currentMixLen;
                    previousY = currentY;
                }
                else if (line.at(1) == 'a') // The age headers
                {
                    in.ignore(maxIgnore, '='); // Kill everything up to the logAge value
                    in >> logAge;
                }
            }
            else // It's an entry line
            {
                int eep;
                double tempMass;

                stringstream in(line);

                vector<double> mags;

                in >> eep >> tempMass;

                // While the stream is still good
                while (in)
                {
                    double mag;
                    in >> mag; // Input a magnitude
                    mags.push_back(mag);
                    in.peek(); // Try the next input to trigger end of stream
                    // This keeps us from getting a duplicate magnitude at the end of the vector.
                }

                if(mags.size() != availableFilters.size())
                {
                    cerr << "Malformed model file:\n\t" << line << endl;
                    exit(1);
                }

                if (!fin.eof())
                {
                    eeps.emplace_back(eep, tempMass, mags);
                }

                haveFirstRecord = true;
            }
        }
    } // EOF

    isochrones.emplace_back(logAge, eeps);
    heliumCurves.emplace_back(previousY, isochrones);
    fehCurves.emplace_back(previousFeH, heliumCurves);

    fin.close();

    ageLimit.first = fehCurves.front().heliumCurves.front().isochrones.front().logAge;
    ageLimit.second = fehCurves.front().heliumCurves.front().isochrones.back().logAge;
}


double GenericMsModel::deriveAgbTipMass (const vector<int> &filters, double newFeH, double newY, double newAge)
{
    isochrone = deriveIsochrone(filters, newFeH, newY, newAge);

    return isochrone.agbTipMass();
}

Isochrone GenericMsModel::deriveIsochrone(const vector<int>& filters, double newFeH, double newY, double newAge) const
{
    if (fehCurves.front().heliumCurves.size() == 1)
        return deriveIsochrone_oneY(filters, newFeH, newAge);
    else
        return deriveIsochrone_manyY(filters, newFeH, newY, newAge);
}


Isochrone GenericMsModel::deriveIsochrone_oneY(const vector<int>& filters, double newFeH, double newAge) const
{
    // Run code comparable to the implementation of deriveAgbTipMass for every mag in every eep, interpolating first in age and then in FeH
    // Check for requested age or [Fe/H] out of bounds
    if ((newAge < ageLimit.first)
     || (newAge > ageLimit.second)
     || (newFeH < fehCurves.front().feh)
     || (newFeH > fehCurves.back().feh))
    {
        // cerr << newAge << " >? " << ageLimit.first << endl;
        // cerr << newAge << " <? " << ageLimit.second << endl;
        // cerr << newFeH << " >? " << fehCurves.front().feh << endl;
        // cerr << newFeH << " <? " << fehCurves.back().feh << endl;
        throw InvalidCluster("Age or FeH out of bounds in GenericMsModel::deriveIsochrone");
    }

    auto fehIter = lower_bound(fehCurves.begin(), fehCurves.end(), newFeH, FehCurve::compareFeh);

    if (fehIter == fehCurves.end())
    {
        fehIter -= 2;
    }
    else if (fehIter != fehCurves.begin())
    {
        fehIter -= 1;
    }

    int iAge;

    {
        auto yIter = fehIter[0].heliumCurves.begin();
        auto ageIter = lower_bound(yIter[0].isochrones.begin(), yIter[0].isochrones.end(), newAge, Isochrone::compareAge);

        if (ageIter == yIter[0].isochrones.end())
        {
            ageIter -= 2;
        }
        else if (ageIter != yIter[0].isochrones.begin())
        {
            ageIter -= 1;
        }

        // Ensure that the ageIter is reasonable
        // Age gridding is identical between Y grid points
        assert(ageIter[0].logAge <  newAge);
        assert(ageIter[1].logAge >= newAge);

        // An age index is more useful, as age indexes should be identical across HeliumCurves
        iAge = ageIter - yIter[0].isochrones.begin();
    }

    vector<Isochrone> interpIso;

    // Interpolate between two ages in two FeHs.
    for (int i = 0; i < 2; ++i)
    {
        // Shortcut iterator
        auto ageIter = fehIter[i].heliumCurves.front().isochrones.begin() + iAge;

        vector<EvolutionaryPoint> interpEeps;

        // Find the minimum and maximum eep values and set a global
        // min and max that is the union of the min and max for each isochrone
        int minEep = std::max(ageIter[0].eeps.front().eep, ageIter[1].eeps.front().eep);
        int maxEep = std::min(ageIter[0].eeps.back().eep, ageIter[1].eeps.back().eep);

        size_t numEeps = maxEep - minEep + 1;     // = the # of eep points that will be in the new isochrone

        // For each isochrone, find the amount by which the min
        // eep # for that isochrone is offset from the global
        // minimum eep #
        array<int, 2> eepOffset = { minEep - ageIter[0].eeps.front().eep
                                  , minEep - ageIter[1].eeps.front().eep };

        for (size_t e = 0; e < numEeps; ++e)
        {
            vector<double> mags;
            mags.resize(FILTS);

            for (auto f : filters)
            {
                mags[f] = linearTransform<>(ageIter[0].logAge
                                          , ageIter[1].logAge
                                          , ageIter[0].eeps.at(e + eepOffset[0]).mags.at(f)
                                          , ageIter[1].eeps.at(e + eepOffset[1]).mags.at(f)
                                          , newAge).val;
            }

            double newMass = linearTransform<>(ageIter[0].logAge
                                             , ageIter[1].logAge
                                             , ageIter[0].eeps.at(e + eepOffset[0]).mass
                                             , ageIter[1].eeps.at(e + eepOffset[1]).mass
                                             , newAge).val;

            // Make sure the EEPs are the same
            assert(ageIter[0].eeps.at(e + eepOffset[0]).eep == ageIter[1].eeps.at(e + eepOffset[1]).eep);

            interpEeps.emplace_back(ageIter[0].eeps.at(e + eepOffset[0]).eep, newMass, mags);
        }


        double interpAge = linearTransform<>(fehIter[0].feh
                                           , fehIter[1].feh
                                           , ageIter[0].logAge
                                           , ageIter[1].logAge
                                           , newFeH).val;

        interpIso.emplace_back(interpAge, interpEeps);
    }

    // Now, interpolate between the two derived isochrones using FeH
    vector<EvolutionaryPoint> interpEeps;
    // Find the minimum and maximum eep values and set a global
    // min and max that is the union of the min and max for each isochrone
    int minEep = std::max(interpIso[0].eeps.front().eep, interpIso[1].eeps.front().eep);
    int maxEep = std::min(interpIso[0].eeps.back().eep, interpIso[1].eeps.back().eep);

    size_t numEeps = maxEep - minEep + 1;     // = the # of eep points that will be in the new isochrone

    // For each isochrone, find the amount by which the min
    // eep # for that isochrone is offset from the global
    // minimum eep #
    array<int, 2> eepOffset = { minEep - interpIso[0].eeps.front().eep
                              , minEep - interpIso[1].eeps.front().eep };

    for (size_t e = 0; e < numEeps; ++e)
    {
        vector<double> mags;
        mags.resize(FILTS);

        for (auto f : filters)
        {
            mags[f] = linearTransform<>(fehIter[0].feh
                                      , fehIter[1].feh
                                      , interpIso.at(0).eeps.at(e + eepOffset[0]).mags.at(f)
                                      , interpIso.at(1).eeps.at(e + eepOffset[1]).mags.at(f)
                                      , newFeH).val;
        }

        double interpMass = linearTransform<>(fehIter[0].feh
                                            , fehIter[1].feh
                                            , interpIso.at(0).eeps.at(e + eepOffset[0]).mass
                                            , interpIso.at(1).eeps.at(e + eepOffset[1]).mass
                                            , newFeH).val;

        // Make sure the EEPs are the same
        assert(interpIso[0].eeps.at(e + eepOffset[0]).eep == interpIso[1].eeps.at(e + eepOffset[1]).eep);

        interpEeps.emplace_back(interpIso.at(0).eeps.at(e + eepOffset[0]).eep, interpMass, mags);
    }

    assert(std::is_sorted(interpEeps.begin(), interpEeps.end()));

    return {interpIso.at(0).logAge, interpEeps};
}


Isochrone GenericMsModel::deriveIsochrone_manyY(const vector<int>& filters, double newFeH, double newY, double newAge) const
{
    // Run code comparable to the implementation of deriveAgbTipMass for every mag in every eep, interpolating first in age and then in FeH
    // Check for requested age or [Fe/H] out of bounds
    if ((newAge < ageLimit.first)
     || (newAge > ageLimit.second)
     || (newFeH < fehCurves.front().feh)
     || (newFeH > fehCurves.back().feh))
    {
        // cerr << newAge << " >? " << ageLimit.first << endl;
        // cerr << newAge << " <? " << ageLimit.second << endl;
        // cerr << newFeH << " >? " << fehCurves.front().feh << endl;
        // cerr << newFeH << " <? " << fehCurves.back().feh << endl;
        throw InvalidCluster("Age or FeH out of bounds in GenericMsModel::deriveIsochrone");
    }

    auto fehIter = lower_bound(fehCurves.begin(), fehCurves.end(), newFeH, FehCurve::compareFeh);

    if (fehIter == fehCurves.end())
    {
        fehIter -= 2;
    }
    else if (fehIter != fehCurves.begin())
    {
        fehIter -= 1;
    }

    int iY;
    int iAge;

    {
        auto yIter = lower_bound(fehIter[0].heliumCurves.begin(), fehIter[0].heliumCurves.end(), newY, HeliumCurve::compareY);

        if (yIter == fehIter[0].heliumCurves.end())
        {
            yIter -= 2;
        }
        else if (yIter != fehIter[0].heliumCurves.begin())
        {
            yIter -= 1;
        }

        // An index is more useful, as Y indexes should be identical across FeHCurves
        iY = yIter - fehIter[0].heliumCurves.begin();

        auto ageIter = lower_bound(yIter[0].isochrones.begin(), yIter[0].isochrones.end(), newAge, Isochrone::compareAge);

        if (ageIter == yIter[0].isochrones.end())
        {
            ageIter -= 2;
        }
        else if (ageIter != yIter[0].isochrones.begin())
        {
            ageIter -= 1;
        }

        // Ensure that the ageIter is reasonable
        // Age gridding is identical between Y grid points
        assert(ageIter[0].logAge <  newAge);
        assert(ageIter[1].logAge >= newAge);

        // An age index is more useful, as age indexes should be identical across HeliumCurves
        iAge = ageIter - yIter[0].isochrones.begin();
    }

    vector<Isochrone> interpIso;

    // Interpolate between two ages in two FeHs.
    for (int i = 0; i < 2; ++i)
    {
        vector<Isochrone> tIso;
        auto yIter = fehIter[i].heliumCurves.begin() + iY;

        for (int y = 0; y < 2; ++y)
        {
            // Shortcut iterator
            auto ageIter = yIter[y].isochrones.begin() + iAge;

            vector<EvolutionaryPoint> interpEeps;

            // Find the minimum and maximum eep values and set a global
            // min and max that is the union of the min and max for each isochrone
            int minEep = std::max(ageIter[0].eeps.front().eep, ageIter[1].eeps.front().eep);
            int maxEep = std::min(ageIter[0].eeps.back().eep, ageIter[1].eeps.back().eep);

            size_t numEeps = maxEep - minEep + 1;     // = the # of eep points that will be in the new isochrone

            // For each isochrone, find the amount by which the min
            // eep # for that isochrone is offset from the global
            // minimum eep #
            array<int, 2> eepOffset = { minEep - ageIter[0].eeps.front().eep
                                      , minEep - ageIter[1].eeps.front().eep };

            for (size_t e = 0; e < numEeps; ++e)
            {
                vector<double> mags;
                mags.resize(FILTS);

                for (auto f : filters)
                {
                    mags[f] = linearTransform<>(ageIter[0].logAge
                                              , ageIter[1].logAge
                                              , ageIter[0].eeps.at(e + eepOffset[0]).mags.at(f)
                                              , ageIter[1].eeps.at(e + eepOffset[1]).mags.at(f)
                                              , newAge).val;
                }

                double newMass = linearTransform<>(ageIter[0].logAge
                                                 , ageIter[1].logAge
                                                 , ageIter[0].eeps.at(e + eepOffset[0]).mass
                                                 , ageIter[1].eeps.at(e + eepOffset[1]).mass
                                                 , newAge).val;

                // Make sure the EEPs are the same
                assert(ageIter[0].eeps.at(e + eepOffset[0]).eep == ageIter[1].eeps.at(e + eepOffset[1]).eep);

                interpEeps.emplace_back(ageIter[0].eeps.at(e + eepOffset[0]).eep, newMass, mags);
            }

            double interpAge = linearTransform<>(yIter[0].y
                                               , yIter[1].y
                                               , ageIter[0].logAge
                                               , ageIter[1].logAge
                                               , newY).val;

            // These won't necessarily be true if extrapolation is allowed (which it
            // is), but should be be true due to the age and FeH checks at the
            // beginning of this function
            assert(ageIter[0].logAge <= interpAge);
            assert(ageIter[1].logAge >= interpAge);

            tIso.emplace_back(interpAge, interpEeps);
        }

        assert(tIso.size() == 2);

        // Now, interpolate between the two derived isochrones using Y
        vector<EvolutionaryPoint> interpEeps;

        // Find the minimum and maximum eep values and set a global
        // min and max that is the union of the min and max for each isochrone
        int minEep = std::max(tIso[0].eeps.front().eep, tIso[1].eeps.front().eep);
        int maxEep = std::min(tIso[0].eeps.back().eep, tIso[1].eeps.back().eep);

        size_t numEeps = maxEep - minEep + 1;     // = the # of eep points that will be in the new isochrone

        // For each isochrone, find the amount by which the min
        // eep # for that isochrone is offset from the global
        // minimum eep #
        array<int, 2> eepOffset = { minEep - tIso[0].eeps.front().eep
                                  , minEep - tIso[1].eeps.front().eep };

        for (size_t e = 0; e < numEeps; ++e)
        {
            vector<double> mags;
            mags.resize(FILTS);

            for (auto f : filters)
            {
                mags[f] = linearTransform<>(yIter[0].y
                                          , yIter[1].y
                                          , tIso.at(0).eeps.at(e + eepOffset[0]).mags.at(f)
                                          , tIso.at(1).eeps.at(e + eepOffset[1]).mags.at(f)
                                          , newY).val;
            }

            double interpMass = linearTransform<>(yIter[0].y
                                                , yIter[1].y
                                                , tIso.at(0).eeps.at(e + eepOffset[0]).mass
                                                , tIso.at(1).eeps.at(e + eepOffset[1]).mass
                                                , newY).val;

            // Make sure the EEPs are the same
            assert(tIso[0].eeps.at(e + eepOffset[0]).eep == tIso[1].eeps.at(e + eepOffset[1]).eep);

            interpEeps.emplace_back(tIso.at(0).eeps.at(e).eep, interpMass, mags);
        }

        double interpAge = linearTransform<>(fehIter[0].feh
                                           , fehIter[1].feh
                                           , tIso[0].logAge
                                           , tIso[1].logAge
                                           , newFeH).val;

        interpIso.emplace_back(interpAge, interpEeps);
    }

    // Now, interpolate between the two derived isochrones using FeH
    vector<EvolutionaryPoint> interpEeps;
    // Find the minimum and maximum eep values and set a global
    // min and max that is the union of the min and max for each isochrone
    int minEep = std::max(interpIso[0].eeps.front().eep, interpIso[1].eeps.front().eep);
    int maxEep = std::min(interpIso[0].eeps.back().eep, interpIso[1].eeps.back().eep);

    size_t numEeps = maxEep - minEep + 1;     // = the # of eep points that will be in the new isochrone

    // For each isochrone, find the amount by which the min
    // eep # for that isochrone is offset from the global
    // minimum eep #
    array<int, 2> eepOffset = { minEep - interpIso[0].eeps.front().eep
                              , minEep - interpIso[1].eeps.front().eep };

    for (size_t e = 0; e < numEeps; ++e)
    {
        vector<double> mags;
        mags.resize(FILTS);

        for (auto f : filters)
        {
            mags[f] = linearTransform<>(fehIter[0].feh
                                      , fehIter[1].feh
                                      , interpIso.at(0).eeps.at(e + eepOffset[0]).mags.at(f)
                                      , interpIso.at(1).eeps.at(e + eepOffset[1]).mags.at(f)
                                      , newFeH).val;
        }

        double interpMass = linearTransform<>(fehIter[0].feh
                                            , fehIter[1].feh
                                            , interpIso.at(0).eeps.at(e + eepOffset[0]).mass
                                            , interpIso.at(1).eeps.at(e + eepOffset[1]).mass
                                            , newFeH).val;

        // Make sure the EEPs are the same
        assert(interpIso[0].eeps.at(e + eepOffset[0]).eep == interpIso[1].eeps.at(e + eepOffset[1]).eep);

        interpEeps.emplace_back(interpIso.at(0).eeps.at(e + eepOffset[0]).eep, interpMass, mags);
    }

    assert(std::is_sorted(interpEeps.begin(), interpEeps.end()));

    return {interpIso.at(0).logAge, interpEeps};
}


vector<double> GenericMsModel::msRgbEvol (const vector<int> &filters, double zamsMass) const
{
    vector<double> mags;
    mags.resize(FILTS);

    auto m = lower_bound(isochrone.eeps.begin(), isochrone.eeps.end(), zamsMass, EvolutionaryPoint::compareMass);

    if (m == isochrone.eeps.end()) {
        m -= 2;
    }
    else if (m != isochrone.eeps.begin()) {
        m -= 1;
    }

    for ( auto f : filters )
    {
        double mag = linearTransform<>(m[0].mass, m[1].mass, m[0].mags.at(f), m[1].mags.at(f), zamsMass).val;

        if (std::fabs(mag) < EPS)
            mags.at(f) = 999.99;
        else
            mags.at(f) = mag;
    }

    return mags;
}


double GenericMsModel::wdPrecLogAge(double thisFeH, double zamsMass, double thisY) const
{
    if (fehCurves.front().heliumCurves.size() == 1)
        return wdPrecLogAge_oneY(thisFeH, zamsMass);
    else
        return wdPrecLogAge_manyY(thisFeH, zamsMass, thisY);
}

double GenericMsModel::wdPrecLogAge_manyY(double thisFeH, double zamsMass, double newY) const
{
    auto fehIter = lower_bound(fehCurves.begin(), fehCurves.end(), thisFeH, FehCurve::compareFeh);

    if (fehIter == fehCurves.end())
    {
        fehIter -= 2;
    }
    else if (fehIter != fehCurves.begin())
    {
        fehIter -= 1;
    }

    int iY;

    {
        auto yIter = lower_bound(fehIter[0].heliumCurves.begin(), fehIter[0].heliumCurves.end(), newY, HeliumCurve::compareY);

        if (yIter == fehIter[0].heliumCurves.end())
        {
            yIter -= 2;
        }
        else if (yIter != fehIter[0].heliumCurves.begin())
        {
            yIter -= 1;
        }

        // An index is more useful, as Y indexes should be identical across FeHCurves
        iY = yIter - fehIter[0].heliumCurves.begin();
    }

    array<double, 2> wdPrecLogAge;

    for (int i = 0; i < 2; ++i)
    {
        auto yIter = fehIter[i].heliumCurves.begin() + iY;

        array<double, 2> tmpPrecLogAge;

        for (int y = 0; y < 2; ++y)
        {
            // The AGBt for the youngest isochrone in the given [Fe/H]
            // This should be the largest AGBt for that [Fe/H]
            if (zamsMass > yIter[y].isochrones.front().agbTipMass())
            {
                tmpPrecLogAge[y] = -2.7 * log10 (zamsMass / yIter[y].isochrones.front().agbTipMass()) + yIter[y].isochrones.front().logAge;
            }
            else
            {
                // Search ages in reverse (since the agbTips decrease as age increases)
                auto ageIter = lower_bound(yIter[y].isochrones.rbegin(), yIter[y].isochrones.rend(), zamsMass, Isochrone::compareAgbTip);

                if (ageIter == yIter[y].isochrones.rend())
                {
                    ageIter -= 2;
                }
                else if (ageIter != yIter[y].isochrones.rbegin())
                {
                    ageIter -= 1;
                }

                // Ensure that we found a reasonable value here
                // This seems backward in the context of the reverse_iterator
                // because the AGBt decreases as logAge increases.

                // This was previously asserted, but it doesn't hold if
                // we're looking for a zamsMass smaller than the minimum
                // (which happens in simCluster)
                // assert(ageIter[0].agbTipMass() <= zamsMass);
                assert(ageIter[1].agbTipMass() > zamsMass);

                tmpPrecLogAge[y] = linearTransform<TransformMethod::Interp>(ageIter[0].agbTipMass()
                                                                          , ageIter[1].agbTipMass()
                                                                          , ageIter[0].logAge
                                                                          , ageIter[1].logAge
                                                                          , zamsMass).val;

                // This seems backwards because of the reverse_iterator
                assert(ageIter[0].logAge >= tmpPrecLogAge[y]);
                assert(ageIter[1].logAge <= tmpPrecLogAge[y]);
            }
        }

        wdPrecLogAge[i] = linearTransform<TransformMethod::Interp>(yIter[0].y
                                                                 , yIter[1].y
                                                                 , tmpPrecLogAge[0]
                                                                 , tmpPrecLogAge[1]
                                                                 , newY).val;

    }

    // Linearly interpolate in FeH
    return linearTransform<TransformMethod::Interp>(fehIter[0].feh, fehIter[1].feh, wdPrecLogAge[0], wdPrecLogAge[1], thisFeH).val;
}

double GenericMsModel::wdPrecLogAge_oneY(double thisFeH, double zamsMass) const
{
    auto fehIter = lower_bound(fehCurves.begin(), fehCurves.end(), thisFeH, FehCurve::compareFeh);

    if (fehIter == fehCurves.end())
    {
        fehIter -= 2;
    }
    else if (fehIter != fehCurves.begin())
    {
        fehIter -= 1;
    }

    array<double, 2> wdPrecLogAge;

    for (int i = 0; i < 2; ++i)
    {
        // The AGBt for the youngest isochrone in the given [Fe/H]
        // This should be the largest AGBt for that [Fe/H]
        if (zamsMass > fehIter[i].heliumCurves.front().isochrones.front().agbTipMass())
        {
            wdPrecLogAge[i] = -2.7 * log10 (zamsMass / fehIter[i].heliumCurves.front().isochrones.front().agbTipMass()) + fehIter[i].heliumCurves.front().isochrones.front().logAge;
        }
        else
        {
            // Search ages in reverse (since the agbTips decrease as age increases)
            auto ageIter = lower_bound(fehIter[i].heliumCurves.front().isochrones.rbegin(), fehIter[i].heliumCurves.front().isochrones.rend(), zamsMass, Isochrone::compareAgbTip);

            if (ageIter == fehIter[i].heliumCurves.front().isochrones.rend())
            {
                ageIter -= 2;
            }
            else if (ageIter != fehIter[i].heliumCurves.front().isochrones.rbegin())
            {
                ageIter -= 1;
            }

            // Ensure that we found a reasonable value here
            // This seems backward in the context of the reverse_iterator
            // because the AGBt decreases as logAge increases.

            // This was previously asserted, but it doesn't hold if
            // we're looking for a zamsMass smaller than the minimum
            // (which happens in simCluster)
            // assert(ageIter[0].agbTipMass() <= zamsMass);
            assert(ageIter[1].agbTipMass() > zamsMass);

            wdPrecLogAge[i] = linearTransform<TransformMethod::Interp>(ageIter[0].agbTipMass()
                                                                     , ageIter[1].agbTipMass()
                                                                     , ageIter[0].logAge
                                                                     , ageIter[1].logAge
                                                                     , zamsMass).val;

            // This seems backwards because of the reverse_iterator
            assert(ageIter[0].logAge >= wdPrecLogAge[i]);
            assert(ageIter[1].logAge <= wdPrecLogAge[i]);
        }
    }

    // Linearly interpolate in FeH
    return linearTransform<TransformMethod::Interp>(fehIter[0].feh, fehIter[1].feh, wdPrecLogAge[0], wdPrecLogAge[1], thisFeH).val;
}
