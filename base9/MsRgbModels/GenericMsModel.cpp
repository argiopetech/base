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
#include "MsRgbModel.hpp"
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
            if ( f == availableFilters[i] )
            {
                indices.push_back(i);

                foundFilter = true;
                break;
            }
        }

        if ( ! foundFilter )
        {
            throw InvalidModelError(f);
        }
    }

    for (size_t i = 0; i < indices.size(); ++i)
    {
        // Prepend the desired filter values to the beginning of the list...
        availableFilters.insert(availableFilters.begin() + i, availableFilters[indices[i] + i]);
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
                        eep.mags.insert(eep.mags.begin() + i, eep.mags[indices[i] + i]);
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


void GenericMsModel::loadModel(string path)
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
        if (!line.empty() && (line[0] != '#'))
        {
            stringstream in(line);

            if (line[0] == '%') // It's a command line
            {
                // Push the previous age
                if (haveFirstRecord)
                {
                    isochrones.emplace_back(logAge, eeps);
                    eeps.clear();
                }

                // A few choices...
                if (line[1] == 'f') // The once-per-file filter list
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
                else if (line[1] == 's') // The section headers
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
                            haveFirstRecord = false;
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
                else if (line[1] == 'a') // The age headers
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


Isochrone* GenericMsModel::deriveIsochrone(double newFeH, double newY, double newAge) const
{
    if (fehCurves.front().heliumCurves.size() == 1)
        return deriveIsochrone_oneY(newFeH, newAge);
    else
        return deriveIsochrone_manyY(newFeH, newY, newAge);
}


Isochrone* GenericMsModel::deriveIsochrone_oneY(double newFeH, double newAge) const
{
    // Run code comparable to the implementation of deriveAgbTipMass for every mag in every eep, interpolating first in age and then in FeH
    // Check for requested age or [Fe/H] out of bounds
    if ((newAge < ageLimit.first)
     || (newAge > ageLimit.second)
     || (newFeH < fehCurves.front().feh)
     || (newFeH > fehCurves.back().feh))
    {
        cerr << newAge << " >? " << ageLimit.first << endl;
        cerr << newAge << " <? " << ageLimit.second << endl;
        cerr << newFeH << " >? " << fehCurves.front().feh << endl;
        cerr << newFeH << " <? " << fehCurves.back().feh << endl;
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

            for (size_t f = 0; f < ageIter[0].eeps[e + eepOffset[0]].mags.size(); ++f)
            {
                mags.push_back(linearTransform<>(ageIter[0].logAge
                                               , ageIter[1].logAge
                                               , ageIter[0].eeps[e + eepOffset[0]].mags[f]
                                               , ageIter[1].eeps[e + eepOffset[1]].mags[f]
                                               , newAge).val);
            }

            double newMass = linearTransform<>(ageIter[0].logAge
                                             , ageIter[1].logAge
                                             , ageIter[0].eeps[e + eepOffset[0]].mass
                                             , ageIter[1].eeps[e + eepOffset[1]].mass
                                             , newAge).val;

            // Make sure the EEPs are the same
            assert(ageIter[0].eeps[e + eepOffset[0]].eep == ageIter[1].eeps[e + eepOffset[1]].eep);

            interpEeps.emplace_back(ageIter[0].eeps[e + eepOffset[0]].eep, newMass, mags);
        }

        interpIso.emplace_back(newAge, interpEeps);
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

        for (size_t f = 0; f < interpIso[0].eeps[e + eepOffset[0]].mags.size(); ++f)
        {
            mags.push_back(linearTransform<>(fehIter[0].feh
                                           , fehIter[1].feh
                                           , interpIso[0].eeps[e + eepOffset[0]].mags[f]
                                           , interpIso[1].eeps[e + eepOffset[1]].mags[f]
                                           , newFeH).val);
        }

        double interpMass = linearTransform<>(fehIter[0].feh
                                            , fehIter[1].feh
                                            , interpIso[0].eeps[e + eepOffset[0]].mass
                                            , interpIso[1].eeps[e + eepOffset[1]].mass
                                            , newFeH).val;

        // Make sure the EEPs are the same
        assert(interpIso[0].eeps[e + eepOffset[0]].eep == interpIso[1].eeps[e + eepOffset[1]].eep);

        interpEeps.emplace_back(interpIso[0].eeps[e + eepOffset[0]].eep, interpMass, mags);
    }

    // This is important, but it takes O(n) time in debug, which is highly un-cool.
    // This should also be checked in marg.cpp when calling calcPost.
    // assert(std::is_sorted(interpEeps.begin(), interpEeps.end()));

    return (new Isochrone(interpIso[0].logAge, interpEeps)) ;
}


Isochrone* GenericMsModel::deriveIsochrone_manyY(double newFeH, double newY, double newAge) const
{
    // Run code comparable to the implementation of deriveAgbTipMass for every mag in every eep, interpolating first in age and then in FeH
    // Check for requested age or [Fe/H] out of bounds
    if ((newAge < ageLimit.first)
     || (newAge > ageLimit.second)
     || (newFeH < fehCurves.front().feh)
     || (newFeH > fehCurves.back().feh))
    {
        cerr << newAge << " >? " << ageLimit.first << endl;
        cerr << newAge << " <? " << ageLimit.second << endl;
        cerr << newFeH << " >? " << fehCurves.front().feh << endl;
        cerr << newFeH << " <? " << fehCurves.back().feh << endl;
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

    bool oneY = false;

    int iY;
    int iAge;

    // Much like wdPrecLogAge, we want to handle the one-Y case properly (with a single function)
    if ( fehIter[0].heliumCurves.size() == 1
      || fehIter[1].heliumCurves.size() == 1 )
    {
        return deriveIsochrone_oneY( newFeH, newAge );

        // Technically unneeded till the two deriveIsochrone functions are combined
        oneY = true;
        iY = 0;
    }
    else
    {
        // This should hold true for all current models which do not trigger the size() == 1 conditional
        assert(fehIter[0].heliumCurves.size() == fehIter[1].heliumCurves.size());

        auto yIter = lower_bound(fehIter->heliumCurves.begin(), fehIter->heliumCurves.end(), newY, HeliumCurve::compareY);

        if (yIter == fehIter->heliumCurves.end())
        {
            yIter -= 2;
        }
        else if (yIter != fehIter->heliumCurves.begin())
        {
            yIter -= 1;
        }

        // An index is more useful, as Y indexes should be identical across FeHCurves
        iY = yIter - fehIter->heliumCurves.begin();
    }

    {
        auto yIter = fehIter->heliumCurves.begin() + iY;
        auto ageIter = lower_bound(yIter->isochrones.begin(), yIter->isochrones.end(), newAge, Isochrone::compareAge);

        if (ageIter == yIter->isochrones.end())
        {
            ageIter -= 2;
        }
        else if (ageIter != yIter->isochrones.begin())
        {
            ageIter -= 1;
        }

        // Ensure that the ageIter is reasonable
        // Age gridding is identical between Y grid points
        assert(ageIter[0].logAge <  newAge);
        assert(ageIter[1].logAge >= newAge);

        // An age index is more useful, as age indexes should be identical across HeliumCurves
        iAge = ageIter - yIter->isochrones.begin();
    }

    vector<Isochrone> interpIso;

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

                for (size_t f = 0; f < ageIter[0].eeps[e + eepOffset[0]].mags.size(); ++f)
                {
                    mags.push_back(linearTransform<>(ageIter[0].logAge
                                                   , ageIter[1].logAge
                                                   , ageIter[0].eeps[e + eepOffset[0]].mags[f]
                                                   , ageIter[1].eeps[e + eepOffset[1]].mags[f]
                                                   , newAge).val);
                }

                double newMass = linearTransform<>(ageIter[0].logAge
                                                 , ageIter[1].logAge
                                                 , ageIter[0].eeps[e + eepOffset[0]].mass
                                                 , ageIter[1].eeps[e + eepOffset[1]].mass
                                                 , newAge).val;

                // Make sure the EEPs are the same
                assert(ageIter[0].eeps[e + eepOffset[0]].eep == ageIter[1].eeps[e + eepOffset[1]].eep);

                interpEeps.emplace_back(ageIter[0].eeps[e + eepOffset[0]].eep, newMass, mags);
            }

            tIso.emplace_back(newAge, interpEeps);
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

            for (size_t f = 0; f < tIso[0].eeps[e + eepOffset[0]].mags.size(); ++f)
            {
                mags.push_back(linearTransform<>(yIter[0].y
                                               , yIter[1].y
                                               , tIso[0].eeps[e + eepOffset[0]].mags[f]
                                               , tIso[1].eeps[e + eepOffset[1]].mags[f]
                                               , newY).val);
            }

            double interpMass = linearTransform<>(yIter[0].y
                                                , yIter[1].y
                                                , tIso[0].eeps[e + eepOffset[0]].mass
                                                , tIso[1].eeps[e + eepOffset[1]].mass
                                                , newY).val;

            // Make sure the EEPs are the same
            assert(tIso[0].eeps[e + eepOffset[0]].eep == tIso[1].eeps[e + eepOffset[1]].eep);

            interpEeps.emplace_back(tIso[0].eeps[e].eep, interpMass, mags);
        }

        interpIso.emplace_back(newAge, interpEeps);
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

        for (size_t f = 0; f < interpIso[0].eeps[e + eepOffset[0]].mags.size(); ++f)
        {
            mags.push_back(linearTransform<>(fehIter[0].feh
                                           , fehIter[1].feh
                                           , interpIso[0].eeps[e + eepOffset[0]].mags[f]
                                           , interpIso[1].eeps[e + eepOffset[1]].mags[f]
                                           , newFeH).val);
        }

        double interpMass = linearTransform<>(fehIter[0].feh
                                            , fehIter[1].feh
                                            , interpIso[0].eeps[e + eepOffset[0]].mass
                                            , interpIso[1].eeps[e + eepOffset[1]].mass
                                            , newFeH).val;

        // Make sure the EEPs are the same
        assert(interpIso[0].eeps[e + eepOffset[0]].eep == interpIso[1].eeps[e + eepOffset[1]].eep);

        interpEeps.emplace_back(interpIso[0].eeps[e + eepOffset[0]].eep, interpMass, mags);
    }

    // assert(std::is_sorted(interpEeps.begin(), interpEeps.end()));
    return (new Isochrone(interpIso[0].logAge, interpEeps));
}


double GenericMsModel::wdPrecLogAge(double thisFeH, double zamsMass, double newY) const
{
    // The guts of the wePrecLogAge function
    // This is implemented internally to avoid namespace pollution and explicit passing of zamsMass
    auto interpY = [=](const vector<HeliumCurve>::const_iterator& yIter)
    {
        double wdPrecLogAge;

        // The AGBt for the youngest isochrone in the given [Fe/H]
        // This should be the largest AGBt for that [Fe/H]
        if (zamsMass > yIter->isochrones.front().agbTipMass())
        {
            wdPrecLogAge = -2.7 * log10 (zamsMass / yIter->isochrones.front().agbTipMass()) + yIter->isochrones.front().logAge;
        }
        else
        {
            // Search ages in reverse (since the agbTips decrease as age increases)
            auto ageIter = lower_bound(yIter->isochrones.rbegin(), yIter->isochrones.rend(), zamsMass, Isochrone::compareAgbTip);

            if (ageIter == yIter->isochrones.rend())
            {
                ageIter -= 2;
            }
            else if (ageIter != yIter->isochrones.rbegin())
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

            wdPrecLogAge = linearTransform<TransformMethod::Interp>(ageIter[0].agbTipMass()
                                                                  , ageIter[1].agbTipMass()
                                                                  , ageIter[0].logAge
                                                                  , ageIter[1].logAge
                                                                  , zamsMass).val;

            // This seems backwards because of the reverse_iterator
            assert(ageIter[0].logAge >= wdPrecLogAge);
            assert(ageIter[1].logAge <= wdPrecLogAge);
        }

        return wdPrecLogAge;
    };

    bool oneY = false;

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

    // If either of the [Fe/H] values has only one Y, we degrade to only running on the first value
    // This is the case in Girardi, Old DSED, and the positive [Fe/H] values of New DSED (at time of writing)
    // Alternatively, we search for a reasonable Y value to interpolate from.
    if ( fehIter[0].heliumCurves.size() == 1
      || fehIter[1].heliumCurves.size() == 1 )
    {
        iY = 0;
        oneY = true;
    }
    else
    {
        // Currently, if neither [Fe/H] has one Y value, they both have the same number of Y values
        assert (fehIter[0].heliumCurves.size() == fehIter[1].heliumCurves.size());

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

        // If we're in the one Y condition from above, this is a simple loop interpolating in [Fe/H]
        // Otherwise, we have to interpolate in Y at both [Fe/H] values
        if (oneY)
        {
            wdPrecLogAge[i] = interpY(yIter);
        }
        else
        {
            array<double, 2> tmpPrecLogAge;

            for (int y = 0; y < 2; ++y)
            {
                tmpPrecLogAge[y] = interpY(yIter + y);
            }

            wdPrecLogAge[i] = linearTransform<TransformMethod::Interp>(yIter[0].y
                                                                     , yIter[1].y
                                                                     , tmpPrecLogAge[0]
                                                                     , tmpPrecLogAge[1]
                                                                     , newY).val;
        }
    }

    // Linearly interpolate in FeH
    return linearTransform<TransformMethod::Interp>(fehIter[0].feh, fehIter[1].feh, wdPrecLogAge[0], wdPrecLogAge[1], thisFeH).val;
}
