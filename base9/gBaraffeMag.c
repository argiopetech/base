#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "evolve.h"
#include "linInterp.h"
#include "gBaraffeMag.h"

// Declared in parent program (simCluster or mcmc, or makeCMD)
extern int    verbose, useFilt[FILTS];
extern double globalMags[FILTS];

static double barAge[N_BAR_AGES];
static double barMass[N_BAR_MASSES];
static double barMag[N_BAR_AGES][N_BAR_MASSES][N_BAR_FILTS];

void loadBaraffe(char *path)
{
  
  int m,a,mStart,filt;
	FILE *pBaraffe;
	char line[1000],tempFile[100]="\0";
  char ind='0';
	double tempAge;
  
	strcat(tempFile, path);
	strcat(tempFile, "COND03_spitzer_0.2-8.2Gyr");
	
  if((pBaraffe = fopen(tempFile,"r")) == NULL) {
    printf("\n\n file %s was not found - exiting\n",tempFile);
    exit(1);
  }
  
	//Skip header lines
  a=-1;
  while(fgets(line,1000,pBaraffe) != NULL){
    sscanf(line," %c ",&ind);
    if(ind == 't'){
      sscanf(line,"%*s %*s %*s %lf ",&tempAge);
      barAge[++a]=log10(tempAge*1e9);
      for(m=0;m<3;m++)  fgets(line, 1000, pBaraffe);
      if(barAge[a]>log10(3.6e9)) mStart = 2;
      else if(barAge[a]>log10(1.4e9)) mStart = 1;
      else mStart = 0;
      for(m=mStart;m<N_BAR_MASSES;m++){
        fscanf(pBaraffe,"%lf %*f %*f %*f ",&barMass[m]);
        for(filt=0;filt<N_BAR_FILTS;filt++) fscanf(pBaraffe,"%lf ",&(barMag[a][m][filt]));
      }
    }
  }
  
  // Kludge to fix missing low mass entries in higher age isochrones
  for(a=23;a<N_BAR_AGES;a++)
    for(filt=0;filt<N_BAR_FILTS;filt++) barMag[a][1][filt] = barMag[a][2][filt] + (barMag[a-1][1][filt] - barMag[a-1][2][filt]);
  for(a=6;a<N_BAR_AGES;a++)
    for(filt=0;filt<N_BAR_FILTS;filt++) barMag[a][0][filt] = barMag[a][1][filt] + (barMag[a-1][0][filt] - barMag[a-1][1][filt]);

  fclose(pBaraffe);
  
}



void getBaraffeMags(double logAge, double mass)
{
	int a, m, filt,i;
	double massMag[2][N_BAR_FILTS];
  
	int    binarySearch (double *searchArray, int size, double searchItem);
  
  a = binarySearch(barAge,N_BAR_AGES,logAge);
  m = binarySearch(barMass,N_BAR_MASSES,mass);
  
  //Interpolate in mass first
  for(i=0;i<2;i++){
    for(filt=0;filt<N_BAR_FILTS;filt++){
      if(useFilt[filt+8]){
        massMag[i][filt] = linInterpExtrap(barMass[m],barMass[m+1], barMag[a+i][m][filt],barMag[a+i][m+1][filt],mass);
      }	
    }
  }
    
  //Interpolate in age
  for(filt=0;filt<N_BAR_FILTS;filt++){
    if(useFilt[filt+8]){
      globalMags[filt+8] = linInterpExtrap(barAge[a],barAge[a+1],massMag[0][filt],massMag[1][filt],logAge);
    }
  }
  for(filt=0;filt<8;filt++) globalMags[filt] = 99.99;
}
