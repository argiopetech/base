#include <iostream>
#include <string>
#include <vector>

#include <cassert>

#include "WdAtmosphereModel.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;

void WdAtmosphereModel::restrictToFilters(const vector<string>& filters)
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
    for (auto &atmos : hCurves)
    {
        for (auto &r : atmos.record)
        {
            for (size_t i = 0; i < indices.size(); ++i)
            {
                // Prepend the desired filter values to the beginning of the list...
                r.mags.insert(r.mags.begin() + i, r.mags.at(indices.at(i) + i));
            }

            // ...then erase the rest.
            r.mags.erase(r.mags.begin() + indices.size(), r.mags.end());
            // And try to reduce the capacity() to what is actually needed
            r.mags.shrink_to_fit();

            // At this point, eep.mags, availableFilters, and filters should all be the same size.
            assert(r.mags.size() <= filters.size());
        }
    }

    for (auto &atmos : heCurves)
    {
        for (auto &r : atmos.record)
        {
            for (size_t i = 0; i < indices.size(); ++i)
            {
                // Prepend the desired filter values to the beginning of the list...
                r.mags.insert(r.mags.begin() + i, r.mags.at(indices.at(i) + i));
            }

            // ...then erase the rest.
            r.mags.erase(r.mags.begin() + indices.size(), r.mags.end());
            // And try to reduce the capacity() to what is actually needed
            r.mags.shrink_to_fit();

            // At this point, eep.mags, availableFilters, and filters should all be the same size.
            assert(r.mags.size() <= filters.size());
        }
    }
}
