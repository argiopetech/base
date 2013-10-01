/*** Proposes new star properties ***/
/*** last update:          ***/

#include <random>

#include "Star.cpp"
#include "samplers.hpp"
#include "densities.hpp"
#include "evolve.hpp"

// void propFieldStar (Star *inputStar)
// {
//     double rand1;

//     rand1 = genrand_res53 ();

//     if (rand1 > (*inputStar).clustStarProposalDens)
//         (*inputStar).isFieldStar = 1;
//     else
//         (*inputStar).isFieldStar = 0;

// }

// void propMass (Star *inputStar)
// {
//     double rand1;

//     rand1 = 2.0 * (*inputStar).UStepSize * ((2.0 * genrand_res53 ()) - 1.0);    // creates uniform random deviate from -2 sigma to +2 sigma
//     (*inputStar).U += rand1;
// }

// void propMassRatio (Star *inputStar)
// {
//     double rand1;

//     if ((*inputStar).status[0] == WD || inputStar->status[0] == BD)
//     {
//         (*inputStar).massRatio = 0.0;   // WD and BD binaries not allowed
//     }
//     else
//     {
//         rand1 = 2.0 * (*inputStar).massRatioStepSize * (2.0 * genrand_res53 () - 1.0);
//         (*inputStar).massRatio += rand1;
//         while ((*inputStar).massRatio < 0.0 || (*inputStar).massRatio > 1.0)
//         {
//             if ((*inputStar).massRatio < 0.0)
//             {
//                 (*inputStar).massRatio = -(*inputStar).massRatio;       //reflect proposed ratio if out of bounds
//             }
//             if ((*inputStar).massRatio > 1.0)
//             {
//                 (*inputStar).massRatio = 2.0 - (*inputStar).massRatio;  //reflect proposed ratio if out of bounds
//             }
//         }
//     }
// }

// void propClustParam (Cluster *clust, int TYPE)
// {
//     if (TYPE == AGE_DURING_WANDER)
//     {
//         (*clust).parameter[AGE] = gen_norm ((*clust).parameter[AGE], 10.0 * AGEPROPSTEPSIZE);
//         return;
//     }
//     else if (TYPE == AGE)
//     {
//         (*clust).parameter[AGE] = (*clust).parameter[AGE] + (*clust).stepSize[AGE] * (genrand_res53 () - 0.5);
//         return;
//     }
//     else
//     {
//         (*clust).parameter[TYPE] = (*clust).parameter[TYPE] + sampleT ((*clust).stepSize[TYPE] * (*clust).stepSize[TYPE]);
//         return;
//     }
// }

double sampleT (std::mt19937 &gen, double var, double nu)
{
    double u = 0.0;
    double y = 0.0;
    double peak = 0.0;

    peak = exp (logTDens (0, 0, var, nu));

    do
    {
        u = std::generate_canonical<double, 53>(gen);

        u = 8.0 * sqrt (var) * (2.0 * u - 1.0);
        y = peak * std::generate_canonical<double, 53>(gen);
        if (y < 1.e-15)
            y = 1.e-15;         /* (TvH) - trap for u = 0; does 53 bit resolution mean 1/(2^53)? */
        y = log (y);
    } while (y > logTDens (u, 0, var, nu));

    return u;
}
