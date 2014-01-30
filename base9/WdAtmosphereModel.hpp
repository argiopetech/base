#ifndef WDATMOS_HPP
#define WDATMOS_HPP

#include <array>
#include <map>
#include <string>
#include <utility>

#include "constants.hpp"

class WdAtmosphereModel
{
  protected:
    struct record
    {
        record(double logTeff, std::array<double, FILTS> mags)
            : logTeff(logTeff), mags(mags)
        {;}

        ~record()
        {;}

        // For use with sort/search functions in <algorithms>
        bool operator<(const struct record &b) const
        {
            return logTeff < b.logTeff;
        }

        double logTeff;

        std::array<double, FILTS> mags;
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

    virtual void loadModel (std::string path, MsFilter filterSet) = 0;
    virtual std::array<double, FILTS> teffToMags  (double wdLogTeff, double wdMass, WdAtmosphere wdType) const = 0;

  protected:
    std::vector<AtmosCurve> hCurves;
    std::vector<AtmosCurve> heCurves;
};
#endif
