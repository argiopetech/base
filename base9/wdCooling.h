#ifndef WDCOOL_H
#define WDCOOL_H

#define MAX_WD_MODEL        150   /* Biggest WD model file (to date) has <~100 entries */

struct wdCarbonCurve {
    double x_carbon;
    int length;
    double logRadius[MAX_WD_MODEL];
    double logAge[MAX_WD_MODEL];
    double logTeff[MAX_WD_MODEL];
};

struct wdCoolingCurve {
    double mass;
    int length;

    double *wdCarbons;
    struct wdCarbonCurve *carbonCurve;
};

void loadWDCool(char *path, int modelSet);
double wdMassToTeffAndRadius(double logAge, double x_carbon, double wdPrecLogAge, double wdMass, double *thisWDLogRadius);
double wdMassToTeffAndRadius_wood(double logAge, double wdPrecLogAge, double wdMass, double *thisWDLogRadius);
double wdMassToTeffAndRadius_montgomery(double logAge, double x_carbon, double wdPrecLogAge, double wdMass, double *thisWDLogRadius);
#endif
