#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cmath>

#include "solve.hpp"

typedef char* FREE_ARG;

const int NR_END = 1;
const double TINY = 1.0e-20;

void solve (double **N, double *yvector, int n)
{

    double d;
    int indx[n];
    int i, j;
    double **a;
    double *b;

    a = matrix (1, n, 1, n);
    b = vector (1, n);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            a[i + 1][j + 1] = N[i][j];
        }
        b[i + 1] = yvector[i];
    }

    ludcmp (a, n, indx, &d);
    lubksb (a, n, indx, b);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            N[i][j] = a[i + 1][j + 1];
        }
        yvector[i] = b[i + 1];
    }
    free_matrix (a, 1, n, 1, n);
    free_vector (b, 1, n);
}


void ludcmp (double **a, int n, int *indx, double *d)
{
    int i, imax = 0, j, k;
    double big, dum, sum, temp;
    double *vv;

    vv = vector (1, n);
    *d = 1.0;
    for (i = 1; i <= n; i++)
    {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if ((temp = fabs (a[i][j])) > big)
                big = temp;
        if (big == 0.0)
            nrerror ("Singular matrix in routine ludcmp");
        vv[i] = 1.0 / big;
    }
    for (j = 1; j <= n; j++)
    {
        for (i = 1; i < j; i++)
        {
            sum = a[i][j];
            for (k = 1; k < i; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i <= n; i++)
        {
            sum = a[i][j];
            for (k = 1; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * fabs (sum)) >= big)
            {
                big = dum;
                imax = i;
            }
        }
        if (j != imax)
        {
            for (k = 1; k <= n; k++)
            {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0)
            a[j][j] = TINY;
        if (j != n)
        {
            dum = 1.0 / (a[j][j]);
            for (i = j + 1; i <= n; i++)
                a[i][j] *= dum;
        }
    }
    free_vector (vv, 1, n);
}

void lubksb (double **a, int n, int *indx, double b[])
{
    int i, ii = 0, ip, j;
    double sum;

    for (i = 1; i <= n; i++)
    {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii)
            for (j = ii; j <= i - 1; j++)
                sum -= a[i][j] * b[j];
        else if (sum)
            ii = i;
        b[i] = sum;
    }
    for (i = n; i >= 1; i--)
    {
        sum = b[i];
        for (j = i + 1; j <= n; j++)
            sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

void nrerror (const char *error_text)
/* Numerical Recipes standard error handler */
{
    fprintf (stderr, "Numerical Recipes run-time error...\n");
    fprintf (stderr, "%s\n", error_text);
    fprintf (stderr, "...now exiting to system...\n");
    exit (1);
}

double *vector (long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    double *v;

    v = (double *) malloc ((size_t) ((nh - nl + 1 + NR_END) * sizeof (double)));
    if (!v)
        nrerror ("allocation failure in vector()");
    return v - nl + NR_END;
}

void free_vector (double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
    free ((FREE_ARG) (v + nl - NR_END));
}

double **matrix (long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double **m;

    /* allocate pointers to rows */
    m = (double **) malloc ((size_t) ((nrow + NR_END) * sizeof (double *)));
    if (!m)
        nrerror ("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (double *) malloc ((size_t) ((nrow * ncol + NR_END) * sizeof (double)));
    if (!m[nrl])
        nrerror ("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

void free_matrix (double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
    free ((FREE_ARG) (m[nrl] + ncl - NR_END));
    free ((FREE_ARG) (m + nrl - NR_END));
}
