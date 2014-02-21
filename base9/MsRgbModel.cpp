#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "binSearch.hpp"
#include "LinearTransform.hpp"
#include "MsRgbModel.hpp"

using std::array;
using std::vector;

array<double, FILTS> MsRgbModel::msRgbEvol (const vector<int> &filters, double zamsMass) const
{
    array<double, FILTS> mags;

    int m2 = binarySearch (isochrone.mass.data(), isochrone.nEntries, zamsMass);

    // auto m = lower_bound(isochrone.mass.begin(), isochrone.mass.end(), zamsMass);;

    // if (m == isochrone.mass.begin()) {
    //     // log << "Poetential mass underflow in MSRGB model" << endl;                                                                                                     
    // }
    // else if (m == isochrone.mass.end()) {
    //     // log << "Mass overflow in MSRGB model" << endl;                                                                                                                 
    //     m -= 2;
    // }
    // else {
    //     m -= 1;
    // }

    for ( auto f : filters )
    {
        if (f < numFilts()) // Do we not like the IFMR models here?
        {
            double mag = linearTransform<>(isochrone.mass.at(m2), isochrone.mass.at(m2 + 1), isochrone.mag.at(m2).at(f), isochrone.mag.at(m2 + 1).at(f), zamsMass).val;

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
