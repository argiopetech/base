#include <random>

#include "samplers.hpp"
#include "densities.hpp"

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
