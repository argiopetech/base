/* leastsquares.h */

#ifdef LEASTSQUARES_H
/* the file has been included already */
#else

#define LEASTSQUARES_H

double leastSquaresBeta (double *x, double *y, int l);
double brokenBeta (double *x, double *y, int l, double *d);
void doubleQuickSort (double *firstx, double *lastx, double *firsty, double *lasty);
void swap (double *p1, double *p2);
void lsqpoly (double *x, double *y, int l, int d, double *c);
void lsqpolyw (double *x, double *y, double *w, int l, int d, double *c);
double powerLaw (double *x, double *y, int l, double *c);

#endif
