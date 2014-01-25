#ifndef WDATMOS_HPP
#define WDATMOS_HPP

#include <array>
#include <map>
#include <string>
#include <utility>

#include "constants.hpp"

class WdCoolingModel
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

  public:
    virtual ~WdCoolingModel() {}


    virtual void loadModel (std::string path, MsFilter filterSet) = 0;
    virtual std::array<double, FILTS> teffToMags  (double wdLogTeff, double wdLogG, WdAtmosphere wdType) const = 0;

  protected:
    std::map<double, std::vector<struct record>> hCurves;
    std::map<double, std::vector<struct record>> heCurves;
};
#endif
