#ifndef WDATMOS_HPP
#define WDATMOS_HPP

#include <map>
#include <string>
#include <utility>

#include "constants.hpp"
#include "StellarModel.hpp"

class WdAtmosphereModel : virtual public StellarModel
{
  protected:
    struct record
    {
        record(double logTeff, double logg, std::vector<double> mags)
            : logTeff(logTeff), logg(logg), mags(mags)
        {;}

        ~record()
        {;}

        static bool compareTeff(const record &a, const double b)
        {
            return a.logTeff < b;
        }

        // For use with sort/search functions in <algorithms>
        bool operator<(const struct record &b) const
        {
            return logTeff < b.logTeff;
        }

        double logTeff;
        double logg;

        std::vector<double> mags;
    };

    struct AtmosCurve
    {
        AtmosCurve(double mass, std::vector<struct record> record)
            : mass(mass), record(record)
        {;}

        // For use with sort/search functions in <algorithms>
        bool operator<(const struct AtmosCurve &b) const
        {
            return mass < b.mass;
        }

        double mass;
        std::vector<struct record> record;
    };

  public:
    virtual ~WdAtmosphereModel() {}

    virtual std::vector<double> teffToMags (double wdLogTeff, double wdMass, WdAtmosphere wdType) const = 0;
    virtual double teffToLogg (double wdLogTeff, double wdMass, WdAtmosphere wdType) const = 0;

    virtual void restrictToFilters(const std::vector<std::string>&) = 0;

    std::vector<std::string> getAvailableFilters() const { return availableFilters; }

  protected:
    std::vector<AtmosCurve> hCurves;
    std::vector<AtmosCurve> heCurves;

    std::vector<std::string> availableFilters;
};
#endif
