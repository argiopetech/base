#if defined( LIN_INTERP_H )
  /* the file has been included already */
#else
#define LIN_INTERP_H

double linInterp(double x1, double x2, double y1, double y2, double xActual);
double linInterpExtrap(double x1, double x2, double y1, double y2, double xActual);

#endif
