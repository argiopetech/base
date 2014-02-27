#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "Cluster.hpp"
#include "LinearTransform.hpp"
#include "MsRgbModel.hpp"

using std::array;
using std::lower_bound;
using std::vector;

array<double, FILTS> MsRgbModel::msRgbEvol (const vector<int> &filters, double zamsMass) const
{
    array<double, FILTS> mags;
    mags.fill(99.999);

    auto m = lower_bound(isochrone.eeps.begin(), isochrone.eeps.end(), zamsMass, EvolutionaryPoint::compareMass);

    if (m == isochrone.eeps.end()) {           
        m -= 2;
    }
    else if (m != isochrone.eeps.begin()) {
        m -= 1;
    }

    for ( auto f : filters )
    {
        if (f < numFilts()) // Do we not like the IFMR models here?
        {
            double mag = linearTransform<>(m[0].mass, m[1].mass, m[0].mags.at(f), m[1].mags.at(f), zamsMass).val;

            if (std::fabs(mag) < EPS)
                mags.at(f) = 999.99;
            else
                mags.at(f) = mag;
        }
        else
        {
            throw std::logic_error("Asked for a filter outside the MSRGB range");
        }
    }

    return mags;
}
