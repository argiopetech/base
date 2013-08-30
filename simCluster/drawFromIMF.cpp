#include <cstdio>
#include <cmath>

#include "mt19937ar.hpp"

/*****************************************************************************************
last update: 28sep05

Draw a star randomly from the Miller-Scalo (1979, ApJS, 41, 513) IMF, using their gaussian
in log(M) form for age = 12 Gyr (c1 = 1.09, c2 = -1.02). Using a mass range of 0.1 to 100
Mo.  (Resetting and recalculating the constants each time this subroutine is called is
inefficient, but more calcuation time will be spent drawing the Gaussian random deviate
(typically done 2+ times per function call) and calculating the final inverse log.)
*****************************************************************************************/
double drawFromIMF ()
{

    double logMass, zamsMass;
    const double mf_sigma = 0.67729, mf_mu = -1.02;     /* sigma = sqrt(1 / 2*c1), mu = c2 */

    double gen_norm (double mean, double std_dev);

    do
    {
        logMass = gen_norm (mf_mu, mf_sigma);
    } while (logMass < -4 /*-0.8238*/  || logMass > 2.0);       /* keep within mass (0.15 + EPS to 100 Mo) limits */

    zamsMass = exp10(logMass);     /* costs about 17% of run time */

    //     log << (" drawFromIMF: zamsMass = %.3f\n", zamsMass);

    return zamsMass;

}
