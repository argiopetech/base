#include <cmath>

#include "poly.hpp"

// Calculates the polynomial c[0] + c[1]*x + c[2]*x^2... c[n]*x^n
double poly (double *c, double x, int n)
{

    double p;
    int j;

    p = c[n];
    for (j = n - 1; j >= 0; j--)
        p = p * x + c[j];

    return p;
}
