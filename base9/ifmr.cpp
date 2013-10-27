#include <iostream>

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "Cluster.hpp"
#include "Star.hpp"
#include "evolve.hpp"
#include "LinearTransform.hpp"

using std::cerr;
using std::endl;

double weidemannIFMR (double zamsMass);
double williamsIFMR (double zamsMass);
double salarisLinearIFMR (double zamsMass);
double salarisPiecewiseIFMR (double zamsMass);
double linearIFMRshift (const Cluster &pCluster, double zamsMass);
double linearIFMR0 (const Cluster &pCluster, double zamsMass);
double linearIFMRage (const Cluster &pCluster, double zamsMass);
double linearIFMRhighShift (const Cluster &pCluster, double zamsMass);
double linearIFMRedit (double zamsMass);
double quadraticIFMRshift (const Cluster &pCluster, double zamsMass);
double quadraticIFMRrotate (const Cluster &pCluster, double zamsMass);
double piecewiseLinearIFMR (const Cluster &pCluster, double zamsMass);

double intlFinalMassReln (const Cluster &pCluster, const Model &evoModels, double zamsMass)
{
    double wdMass = 0.0;

    if (evoModels.IFMR == WEIDEMANN)
        wdMass = weidemannIFMR (zamsMass);
    else if (evoModels.IFMR == WILLIAMS)
        wdMass = williamsIFMR (zamsMass);
    else if (evoModels.IFMR == SALARISLIN)
        wdMass = salarisLinearIFMR (zamsMass);
    else if (evoModels.IFMR == SALARISPW)
        wdMass = salarisPiecewiseIFMR (zamsMass);
    else if (evoModels.IFMR == LINEAR)
        wdMass = linearIFMRshift (pCluster, zamsMass);
    else if (evoModels.IFMR == 5)
        wdMass = linearIFMR0 (pCluster, zamsMass);
    else if (evoModels.IFMR == 6)
        wdMass = linearIFMRage (pCluster, zamsMass);
    else if (evoModels.IFMR == 7)
        wdMass = linearIFMRhighShift (pCluster, zamsMass);
    else if (evoModels.IFMR == 8)
        wdMass = linearIFMRedit (zamsMass);
    else if (evoModels.IFMR == 9)
        wdMass = quadraticIFMRshift (pCluster, zamsMass);
    else if (evoModels.IFMR == 10)
        wdMass = quadraticIFMRrotate (pCluster, zamsMass);
    else if (evoModels.IFMR == 11)
        wdMass = piecewiseLinearIFMR (pCluster, zamsMass);
    else
        cerr << "ERROR: Undefined IFMR" << endl;
    return wdMass;
}

double linearIFMRshift (const Cluster &pCluster, double zamsMass)
{
    constexpr double shiftMass = 3.0;
    double wdMass = pCluster.ifmrIntercept + pCluster.ifmrSlope * (zamsMass - shiftMass);

    return wdMass;
}

double linearIFMRage (const Cluster &pCluster, double zamsMass)
{
    double wdMass = pCluster.ifmrIntercept + pCluster.ifmrSlope * (zamsMass - 27.0 + 2.6 * pCluster.age);

    return wdMass;
}

double linearIFMR0 (const Cluster &pCluster, double zamsMass)
{
    double wdMass = pCluster.ifmrIntercept + pCluster.ifmrSlope * zamsMass;

    return wdMass;
}

double linearIFMRhighShift (const Cluster &pCluster, double zamsMass)
{
    constexpr double shiftMass = 5.0;
    double wdMass = pCluster.ifmrIntercept + pCluster.ifmrSlope * (zamsMass - shiftMass);

    return wdMass;
}

double linearIFMRedit (double zamsMass)
{
    constexpr double shiftMass = 2.0;
    double wdMass = 0.985 + 0.13 * (zamsMass - shiftMass);

    return wdMass;
}

double williamsIFMR (double zamsMass)
{
    /***************************************************************************************
        last update: 02aug10
        Linear IFMR from Williams, Bolte, & Koester (2009) ApJ 693:355

        M_final = 0.339 +- 0.015 + (0.129 +- 0.004)*M_init
    ***************************************************************************************/

    double wdMass = 0.0;

    wdMass = 0.339 + 0.129 * zamsMass;
    return wdMass;
}


double salarisLinearIFMR (double zamsMass)
{
    /***************************************************************************************
        last update: 02aug10
        Linear IFMR from Salaris, et al. (2009) ApJ 692:1013

        M_final = 0.466 + 0.084 * M_init
    ***************************************************************************************/

    double wdMass = 0.0;

    wdMass = 0.466 + 0.084 * zamsMass;
    return wdMass;
}


double salarisPiecewiseIFMR (double zamsMass)
{
    /***************************************************************************************
        last update: 02aug10
        Linear IFMR from Salaris, et al. (2009) ApJ 692:1013

        M_final = 0.331 + 0.134 * M_init   1.7_sun < M_init < 4 M_sun
                  0.679 + 0.047 * M_init   4 M_sun < M_init
    ***************************************************************************************/

    double wdMass = 0.0;

    if (zamsMass > 4)
        wdMass = 0.679 + 0.047 * zamsMass;
    else
        wdMass = 0.331 + 0.134 * zamsMass;
    return wdMass;
}

double quadraticIFMRshift (const Cluster &pCluster, double zamsMass)
{
    constexpr double shiftMass = 3.0;
    double wdMass = pCluster.ifmrIntercept + pCluster.ifmrSlope * (zamsMass - shiftMass) + pCluster.ifmrQuadCoef * (zamsMass - shiftMass) * (zamsMass - shiftMass);

    return wdMass;
}

double quadraticIFMRrotate (const Cluster &pCluster, double zamsMass)
{
    constexpr double shiftMass = 3.0;
    constexpr double massLower = 0.15;
    double massUpper = pCluster.M_wd_up;

    double angle = atan (pCluster.ifmrSlope);

    double xLo = (cos (angle) + pCluster.ifmrSlope * sin (angle)) * (massLower - shiftMass);
    double xUp = (cos (angle) + pCluster.ifmrSlope * sin (angle)) * (massUpper - shiftMass);

    /* constants to simplify expressions */
    double cLo = 1.0 / tan (angle) * (zamsMass - shiftMass) - 1.0 / sin (angle) * xLo;
    double cUp = 1.0 / tan (angle) * (zamsMass - shiftMass) - 1.0 / sin (angle) * xUp;
    double A = cLo + cUp - cos (angle) / (sin (angle) * sin (angle) * pCluster.ifmrQuadCoef);

    double wdMass = pCluster.ifmrIntercept - 0.5 * A + sqrt (0.25 * A * A - cLo * cUp - 1.0 / (sin (angle) * pCluster.ifmrQuadCoef) * (zamsMass - shiftMass));

    return wdMass;
}

double piecewiseLinearIFMR (const Cluster &pCluster, double zamsMass)
{
    double shiftMass = 3.0;
    double breakpointMass = 4.0;
    double wdMass;

    if (zamsMass < breakpointMass)
    {
        wdMass = pCluster.ifmrIntercept + pCluster.ifmrSlope * (zamsMass - shiftMass);
    }
    else
    {
        wdMass = pCluster.ifmrIntercept + pCluster.ifmrSlope * (breakpointMass - shiftMass) + pCluster.ifmrQuadCoef * (zamsMass - breakpointMass);
    }

    return wdMass;
}

double weidemannIFMR (double zamsMass)
/***************************************************************************************
last update: 18sep05

Interpolate (1-D) within the Weidemann (Weidemann, V. 2000, A&A, 363, 647) initial-final
mass relation.  This is consistent with the relation presented by Claver, C.F., Liebert,
J, Bergeron, P., & Koester, D. 2001, ApJ, 563, 987 (see figure 11), except the latter
authors' relation has Mf ~ 0.1 Mo higher after Mi = 3 Mo than the former relation.
This can be simulated later by boosting final masses from >= 3 Mo ZAMS stars by 0.1 Mo.

The relation given by Weidemann is

ZAMS  1.0  2.0  2.5  3.0  4.0  5.0  6.0  7.0
WD    0.55 0.60 0.63 0.68 0.80 0.88 0.95 1.02

Extend the Weidemann et al. relation artificially to 9 Mo, so that M_wd_up can be as
large as 9.0.
***************************************************************************************/
{
    constexpr std::array<double, 10> initMassW = {{1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}};
    constexpr std::array<double, 10> finalMassW = {{0.55, 0.60, 0.63, 0.68, 0.80, 0.88, 0.95, 1.02, 1.2, 1.4}};

    if (zamsMass <= initMassW[0])
        return finalMassW[0];   // if out of init-final calib range, use
    else if (zamsMass >= initMassW[9])
        return finalMassW[9];   // limiting values
    else
    {
        int lo = 1;
        int hi = 10;
        int mid;

        // binary search on zams_mass
        while (1)
        {
            mid = ((lo + hi) >> 1);

            if (initMassW[mid - 1] <= zamsMass && zamsMass <= initMassW[mid])
            {
                auto wdMass = linearTransform<TransformMethod::Interp>(initMassW[mid - 1], initMassW[mid], finalMassW[mid - 1], finalMassW[mid], zamsMass);
                return wdMass.val;
            }

            if (lo >= hi)
            {
                cerr << "ERROR: BINARY SEARCH FAILURE" << endl;
                break;
            }

            if (zamsMass > initMassW[mid])
                lo = mid + 1;
            if (zamsMass < initMassW[mid])
                hi = mid - 1;
        }
    }
    return 0.0;
}
