//==================================================== file = gennorm.c =====
//=  Subroutine to generate nomrally distributed random variables           =
//=                                                                         =
//=  Uses the Box-Muller method and only generates one of two paired normal =
//=  random variables.                                                      =
//===========================================================================
//=-------------------------------------------------------------------------=
//=  Build: bcc32 gennorm.c                                                 =
//=-------------------------------------------------------------------------=
//=  Author: Kenneth J. Christensen                                         =
//=          University of South Florida                                    =
//=          WWW: http://www.csee.usf.edu/~christen                         =
//=          Email: christen@csee.usf.edu                                   =
//=-------------------------------------------------------------------------=
//=  History: KJC (06/06/02) - Genesis                                      =
//=           KJC (05/20/03) - Added Jain's RNG for finer granularity       =
//=           TvH (09/28/05) - use norm (now gen_norm) subroutine only, use =
//=                            genrand_res53() instead of Jain's RNG,       =
//=                            rand_val().  note: need to run               =
//=                            "init_genrand(seed);" before invoking this   =
//=                            subroutine                                   =
//===========================================================================

#include <math.h>               // Needed for sqrt() and log()
#include "mt19937ar.h"

#define PI         3.14159265   // The value of pi

//----- Function prototypes -------------------------------------------------
double gen_norm (double mean, double std_dev);  // Returns a normal rv

//===========================================================================
//=  Function to generate normally distributed random variable using the    =
//=  Box-Muller method                                                      =
//=    - Input: mean and standard deviation                                 =
//=    - Output: Returns with normally distributed random variable          =
//===========================================================================
double gen_norm (double mean, double std_dev)
{
    double u, r, theta;         // Variables for Box-Muller method
    double x;                   // Normal(0, 1) rv
    double norm_rv;             // The adjusted normal rv

    // Generate u
    u = 0.0;
    while (u == 0.0)
        u = genrand_res53 ();

    // Compute r
    r = sqrt (-2.0 * log (u));

    // Generate theta
    theta = 0.0;
    while (theta == 0.0)
        theta = 2.0 * PI * genrand_res53 ();

    // Generate x value
    x = r * cos (theta);

    // Adjust x value for specified mean and variance
    norm_rv = (x * std_dev) + mean;

    // Return the normally distributed RV value
    return (norm_rv);
}
