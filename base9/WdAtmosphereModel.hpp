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
        record(double logG, std::array<double, FILTS> mags)
            : logG(logG), mags(mags)
        {;}

        ~record()
        {;}

        // For use with sort/search functions in <algorithms>
        bool operator<(const struct record &b) const
        {
            return logG < b.logG;
        }

        double logG;

        std::array<double, FILTS> mags;
    };

    struct AtmosCurve
    {
        AtmosCurve(double teff, std::vector<struct record> record)
            : logTeff(teff), record(record)
        {;}

        // For use with sort/search functions in <algorithms>
        bool operator<(const struct AtmosCurve &b) const
        {
            return logTeff < b.logTeff;
        }

        double logTeff;
        std::vector<struct record> record;
    };

  public:
    virtual ~WdAtmosphereModel() {}

    virtual void loadModel (std::string path, MsFilter filterSet) = 0;
    virtual std::array<double, FILTS> teffToMags  (double wdLogTeff, double wdLogG, WdAtmosphere wdType) const = 0;

  protected:
    std::vector<AtmosCurve> hCurves;
    std::vector<AtmosCurve> heCurves;
};
#endif
