#include <iostream>

#include <cmath>

#include "Cluster.hpp"
#include "Star.hpp"
#include "LinearTransform.hpp"

using std::cerr;
using std::endl;

double weidemannIFMR (double);
double williamsIFMR (double);
double salarisLinearIFMR (double);
double salarisPiecewiseIFMR (double);
double linearIFMRedit (double);
double linearIFMRshift (const Cluster&, double);
double linearIFMR0 (const Cluster&, double);
double linearIFMRage (const Cluster&, double);
double linearIFMRhighShift (const Cluster&, double);
double quadraticIFMRshift (const Cluster&, double);
double quadraticIFMRrotate (const Cluster&, double);
double piecewiseLinearIFMR (const Cluster&, double);

double intlFinalMassReln (const Cluster &clust, const Model &evoModels, double zamsMass)
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
        wdMass = linearIFMRshift (clust, zamsMass);
    else if (evoModels.IFMR == 5)
        wdMass = linearIFMR0 (clust, zamsMass);
    else if (evoModels.IFMR == 6)
        wdMass = linearIFMRage (clust, zamsMass);
    else if (evoModels.IFMR == 7)
        wdMass = linearIFMRhighShift (clust, zamsMass);
    else if (evoModels.IFMR == 8)
        wdMass = linearIFMRedit (zamsMass);
    else if (evoModels.IFMR == 9)
        wdMass = quadraticIFMRshift (clust, zamsMass);
    else if (evoModels.IFMR == 10)
        wdMass = quadraticIFMRrotate (clust, zamsMass);
    else if (evoModels.IFMR == 11)
        wdMass = piecewiseLinearIFMR (clust, zamsMass);
    else
        cerr << "ERROR: Undefined IFMR" << endl;
    return wdMass;
}

double linearIFMRshift (const Cluster &clust, double zamsMass)
{
    constexpr double shiftMass = 3.0;
    double wdMass = clust.ifmrIntercept + clust.ifmrSlope * (zamsMass - shiftMass);

    return wdMass;
}

double linearIFMRage (const Cluster &clust, double zamsMass)
{
    double wdMass = clust.ifmrIntercept + clust.ifmrSlope * (zamsMass - 27.0 + 2.6 * clust.age);

    return wdMass;
}

double linearIFMR0 (const Cluster &clust, double zamsMass)
{
    double wdMass = clust.ifmrIntercept + clust.ifmrSlope * zamsMass;

    return wdMass;
}

double linearIFMRhighShift (const Cluster &clust, double zamsMass)
{
    constexpr double shiftMass = 5.0;
    double wdMass = clust.ifmrIntercept + clust.ifmrSlope * (zamsMass - shiftMass);

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

double quadraticIFMRshift (const Cluster &clust, double zamsMass)
{
    constexpr double shiftMass = 3.0;
    double wdMass = clust.ifmrIntercept + clust.ifmrSlope * (zamsMass - shiftMass) + clust.ifmrQuadCoef * (zamsMass - shiftMass) * (zamsMass - shiftMass);

    return wdMass;
}

double quadraticIFMRrotate (const Cluster &clust, double zamsMass)
{
    constexpr double shiftMass = 3.0;
    constexpr double massLower = 0.15;
    double massUpper = clust.getM_wd_up();

    double angle = atan (clust.ifmrSlope);

    double xLo = (cos (angle) + clust.ifmrSlope * sin (angle)) * (massLower - shiftMass);
    double xUp = (cos (angle) + clust.ifmrSlope * sin (angle)) * (massUpper - shiftMass);

    /* constants to simplify expressions */
    double cLo = 1.0 / tan (angle) * (zamsMass - shiftMass) - 1.0 / sin (angle) * xLo;
    double cUp = 1.0 / tan (angle) * (zamsMass - shiftMass) - 1.0 / sin (angle) * xUp;
    double A = cLo + cUp - cos (angle) / (sin (angle) * sin (angle) * clust.ifmrQuadCoef);

    double wdMass = clust.ifmrIntercept - 0.5 * A + sqrt (0.25 * A * A - cLo * cUp - 1.0 / (sin (angle) * clust.ifmrQuadCoef) * (zamsMass - shiftMass));

    return wdMass;
}

double piecewiseLinearIFMR (const Cluster &clust, double zamsMass)
{
    double shiftMass = 3.0;
    double breakpointMass = 4.0;
    double wdMass;

    if (zamsMass < breakpointMass)
    {
        wdMass = clust.ifmrIntercept + clust.ifmrSlope * (zamsMass - shiftMass);
    }
    else
    {
        wdMass = clust.ifmrIntercept + clust.ifmrSlope * (breakpointMass - shiftMass) + clust.ifmrQuadCoef * (zamsMass - breakpointMass);
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
        for (int i = 0; i < 10; ++i)
        {
            if ((initMassW[i - 1] <= zamsMass) && (zamsMass <= initMassW[i]))
            {
                return linearTransform<TransformMethod::Interp>(initMassW[i - 1], initMassW[i], finalMassW[i - 1], finalMassW[i], zamsMass).val;
            }
        }
    }

    return 0.0;
}
