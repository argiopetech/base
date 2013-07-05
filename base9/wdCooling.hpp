#ifndef WDCOOL_H
#define WDCOOL_H

#include <string>
#include <vector>

struct record
{
    record(double logRadius, double logAge, double logTeff)
        : logRadius(logRadius), logAge(logAge), logTeff(logTeff)
    {;}

    ~record()
    {;}

    // These are for use with search/sort functions in e.g., <algorithms>
    static bool compareAge(const record &a, const record &b)
    {
        return a.logAge < b.logAge;
    }

    static bool compareRadius(const record &a, const record &b)
    {
        return a.logRadius < b.logRadius;
    }

    static bool compareTeff(const record &a, const record &b)
    {
        return a.logTeff < b.logTeff;
    }

    double logRadius;
    double logAge;
    double logTeff;
};

struct wdCarbonCurve
{
    wdCarbonCurve(double carbon)
        : carbon(carbon), records()
    {;}

    ~wdCarbonCurve()
    {;}

    bool operator<(const struct wdCarbonCurve &b) const
    {
        return carbon < b.carbon;
    }

    const double carbon;

    std::vector<record> records;
};

struct wdCoolingCurve
{
    wdCoolingCurve(double mass)
        : mass(mass), carbonCurves()
    {;}

    ~wdCoolingCurve()
    {;}

    bool operator<(const struct wdCoolingCurve &b) const
    {
        return mass < b.mass;
    }

    const double mass;

    std::vector<struct wdCarbonCurve> carbonCurves;
};

void loadWDCool (std::string path, int modelSet);
double wdMassToTeffAndRadius (double logAge, double x_carbon, double wdPrecLogAge, double wdMass, double &thisWDLogRadius);
#endif
