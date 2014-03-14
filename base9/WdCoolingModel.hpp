#ifndef WDCOOLING_HPP
#define WDCOOLING_HPP

#include <string>
#include <utility>
#include <vector>

#include "StellarModel.hpp"

class WdCoolingModel : public StellarModel
{
  protected:
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

        virtual ~wdCoolingCurve()
        {;}

        bool operator<(const struct wdCoolingCurve &b) const
        {
            return mass < b.mass;
        }

        const double mass;

        std::vector<struct wdCarbonCurve> carbonCurves;
    };

  public:
    virtual ~WdCoolingModel() {}

    virtual std::pair<double, double> wdMassToTeffAndRadius (double logAge, double x_carbon, double wdPrecLogAge, double wdMass) const = 0;

    // WD Cooling models are FilterSet agnostic
    virtual bool isSupported(FilterSetName) { return true; }

  protected:
    std::vector<struct wdCoolingCurve> wdCurves;
};
#endif
