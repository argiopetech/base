#include <stdio.h>
#include <math.h>
#include "linInterp.hpp"

const double EPS = 1e-6;

/*****************************************************************************************
last update: 15sep05

For a given grid with values x1, y1 and x2, y2 with a value of interest, xActual, linearly
interpolate to determine yActual.
*****************************************************************************************/
double linInterp (double x1, double x2, double y1, double y2, double xActual)
{
    double position;

    if (fabs (x2 - x1) < EPS)
    {
        return y1;                      /* which should equal y2 */
    }
    else if (((xActual < x1) && (xActual < x2)) || ((xActual > x1) && (xActual > x2)))
    {
        return 0.0;                     /* use for unknown value */
    }
    else
    {
        position = (xActual - x1) / (x2 - x1);  /* now safe from div by 0 */
        return y1 + position * (y2 - y1);       /* this is yActual */
    }
}

/*****************************************************************************************
Will extrapolate to a point out of range.  This is essentially the same subroutine as
lin_interp(), except that it does not trap for extrapolation, and therefore the calcualted
value can be either linear interpolation or linear extrapolation.
*****************************************************************************************/
double linInterpExtrap (double x1, double x2, double y1, double y2, double xActual)
{
    double position;

    if (fabs (x2 - x1) < EPS)
    {
        return y1;                      /* which should equal y2 */
    }
    else
    {
        position = (xActual - x1) / (x2 - x1);  /* now safe from div by 0 */
        return y1 + position * (y2 - y1);       /* this is yActual */
    }

}
