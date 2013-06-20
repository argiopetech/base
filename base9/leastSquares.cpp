/*** last update:   29aug10  ***/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "leastSquares.h"
#include "poly.h"
#include "solve.h"

// Calculates beta using least squares fit (in order to eliminate correlation between ages and WD masses)
double leastSquaresBeta(double *x, double *y, int l)
{
   double xSum, ySum, num, den;
   int i;
   xSum = 0.0;
   ySum = 0.0;
   num = 0.0;
   den = 0.0;
  
   for(i=0; i < l; i++)
   {
      xSum += x[i];
      ySum += y[i];
      num += (x[i] * y[i]);
      den += (x[i] * x[i]);
   }

   num *= (double) l;
   den *= (double) l;
 
   num -= xSum * ySum;
   den -= xSum * xSum;

   return (double) num / den;
}

// Calculates the best fit singly-broken line to
// a set of data.  Returns the break point and 
// stores the the two slopes in b (which should be 
// a 2-D array
double brokenBeta(double *x, double *y, int l, double *d){
  //Not complete--don't use
  //Need to write output
  int i,bin,lo, n,j;
  double lowchi,ay;  
  double chi[l], xy[l],B1[l],B2[l];

  for(i=0;i<l;i++) chi[i] = xy[i] = B1[i] = B2[i] = 0.0;

  //Divide into 100 bins
  bin = (int) (l/100.0);

  // Sort the arrays in ascending order in x
  doubleQuickSort(x,&(x[l-1]),y,&(y[l-1]));

  j=-1;
  lowchi = HUGE_VAL;

  for(i=l/10.0;i<9.0*l/10.0;i++){

    // Calculate the slopes of the correlations 
    // before and after the pivot point
    B1[i] = leastSquaresBeta(&(x[0]),&(y[0]),i+1);
    B2[i] = leastSquaresBeta(&(x[i+1]),&(y[i+1]),l-(i+1));

    // Find the average value of y at the pivot point
    if(i < bin/2.0) lo = 0;
    else if(i >= l - bin/2.0) lo = l-bin;
    else lo=i-bin/2.0;
    ay = 0.0;
    for(n=lo;n<lo + bin;n++) ay += y[n];
    ay /= bin;

    // Calculate the expected value of y at point x
    // and compare with the actual value to get chi^2
    for(n=0;n<i+1;n++){
      xy[n] = B1[i]*(x[n]-x[i])+ay;
      chi[i] += (xy[n] - y[n])*(xy[n] - y[n]);
    }
    for(n=i+1;n<l;n++){
      xy[n] = B2[i]*(x[n]-x[i])+ay;
      chi[i] += (xy[n] - y[n])*(xy[n] - y[n]);
    }

    //printf("%d %f %f %f %f %f %f\n",i,x[i],y[i],B1[i],B2[i],ay,chi[i]);

    if(chi[i] < lowchi){
      lowchi = chi[i];
      j=i;
    }
  }

  d[0] = B1[j];
  d[1] = B2[j];
  return x[j];

  return 0.0;
}

/*
double findChiSq(double *x, double *y, int l,double a){

  int i;

  i = binarySearch(x,l,a);


}
*/

void swap(double *p1, double *p2){
  double temp;

  temp = *p1;
  *p1 = *p2;
  *p2 = temp;
}

void doubleQuickSort(double *firstx, double *lastx, double *firsty, double *lasty){

  double *leftx = firstx;
  double *rightx = lastx;
  double *lefty = firsty;
  double *righty = lasty;
  double pivot_value;

  ptrdiff_t p_diff;

  p_diff = (rightx - leftx) >> 1;
  pivot_value = *(leftx + p_diff);

  while (leftx <= rightx){
    while (*leftx < pivot_value){
      leftx++;
      lefty++;
    }
    while (*rightx > pivot_value){
      rightx--;
      righty--;
    }
    if(leftx <= rightx){
      swap(leftx, rightx);
      swap(lefty,righty);
      leftx++;
      lefty++;
      rightx--;
      righty--;
    }
  }

  if(firstx < rightx)
    doubleQuickSort(firstx,rightx,firsty,righty);
  if(leftx < lastx)
    doubleQuickSort(leftx,lastx,lefty,lasty);
}

// Method for determining the best fit polynomial of degree d 
// to data x and y (of length l), given weights w and stores the
// results in c (which must have at least d+1 elements)
void lsqpolyw(double *x, double *y, double *w, int l, int d, double *c){

  int i,j,k;
  double **N;//[d+1][d+1];
  double yvector[d+1];

  if((N = (double **) calloc(d+1, sizeof(double*))) == NULL)
    perror("MEMORY ALLOCATION ERROR \n");
  for(j = 0; j < d+1; j++) {
    if((N[j] = (double *) calloc(d+1, sizeof(double))) == NULL)
      perror("MEMORY ALLOCATJON ERROR \n");
  }

  for(i=0;i<d+1;i++){
    for(j=0;j<d+1;j++){
      N[i][j]=0.0;
      for(k=0;k<l;k++) N[i][j] += w[k]*pow(x[k],i+j);
    }
    yvector[i] = 0.0;
    for(k=0;k<l;k++) yvector[i] += w[k]*y[k]*pow(x[k],i);
  }

  solve(N, yvector, d+1);
  for(i=0;i<d+1;i++) c[i] = yvector[i];
}

// Method for determining the best fit polynomial of degree d 
// to data x and y (of length l), with all weights = 1 and stores the
// results in c
void lsqpoly(double *x, double *y, int l, int d, double *c){

  int i,j,k;
  double **N;
  double yvector[d+1];
  
  if((N = (double **) calloc(d+1, sizeof(double*))) == NULL)
    perror("MEMORY ALLOCATION ERROR \n");
  for(j = 0; j < d+1; j++) {
    if((N[j] = (double *) calloc(d+1, sizeof(double))) == NULL)
      perror("MEMORY ALLOCATJON ERROR \n");
  }
  for(i=0;i<d+1;i++){
    for(j=0;j<d+1;j++){
      for(k=0;k<l;k++){
        N[i][j] += pow(x[k],i+j);
      }
    }
    yvector[i] = 0.0;
    for(k=0;k<l;k++) yvector[i] += y[k]*pow(x[k],i);
  }
  if(N[1][1]< 1e-9){
    for(i=0;i<d+1;i++) c[i] =0.0;
    return;
  }
  solve(N, yvector, d+1);
  for(i=0;i<d+1;i++) c[i] = yvector[i];
}

// Calculates the best fit power law between x and y
// such that y = a[k][1]*x^p[k] = c*x^p[k]
double powerLaw(double *x,double *y, int l, double *c){

  int i,j,k=0;
  // double lo=0.01, hi=20,tol =0.01,minchisq;
  double lo=1, hi=20, tol=0.01, minchisq; // (NS) setting the minimum power to 1
  double p[5], a[5][2],chisq[5], xx[l];

  while(hi-lo>tol){

    // Calculate the best fit slope and intercept for (x^p[i],y)
    // for 5 different values of p
    for(i=0;i<5;i++){
      p[i] = lo + i*(hi-lo)/4.;
      chisq[i] = a[i][0] = a[i][1] = 0.0;
      for(j=0;j<l;j++) xx[j] = pow(x[j],p[i]);
      lsqpoly(xx,y,l,1,a[i]);
      for(j=0;j<l;j++) chisq[i]+=pow(y[j]-poly(a[i],xx[j],1),2);
    }

    // Find the value of p with the lowest chisq
    minchisq = HUGE_VAL;
    k=-1;
    for(i=0;i<5;i++){
      if(chisq[i]<minchisq){
        minchisq=chisq[i];
        k=i;
      }
    }

    // Set new hi and low values to calculate another 5 p values
    if(k==1) hi=p[2];
    else if (k==5) lo=p[3];
    else{
      lo=p[k-1];
      hi=p[k+1];
    }
  }

  *c = a[k][1];
  return p[k];


}
