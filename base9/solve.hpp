#ifndef SOLVE_H
#define SOLVE_H

void nrerror (const char *error_text);
double *vector (long nl, long nh);
void free_vector (double *v, long nl, long nh);
void ludcmp (double **a, int n, int *indx, double *d);
void lubksb (double **a, int n, int *indx, double b[]);
void solve (double **a, double *b, int n);
double **matrix (long nrl, long nrh, long ncl, long nch);
void free_matrix (double **m, long nrl, long nrh, long ncl, long nch);

#endif
