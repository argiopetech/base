#ifndef MSRGBEVOL_HPP
#define MSRGBEVOL_HPP

#include <string>
#include <utility>
#include <vector>

#include "constants.hpp"
#include "Isochrone.hpp"
#include "StellarModel.hpp"

class MsRgbModel : virtual public StellarModel
{
  protected:
    struct HeliumCurve
    {
        HeliumCurve(double y, std::vector<Isochrone> isochrones)
            : y(y), isochrones(isochrones)
        {;}

        ~HeliumCurve()
        {;}

        static bool compareY(const HeliumCurve &a, const double b)
        {
            return a.y < b;
        }

        bool operator<(const struct HeliumCurve &b) const
        {
            return y < b.y;
        }

        double y;
        std::vector<Isochrone> isochrones;
    };

    struct FehCurve
    {
        FehCurve(double feh, std::vector<HeliumCurve> heliumCurves)
            : feh(feh), heliumCurves(heliumCurves)
        {;}

        ~FehCurve()
        {;}

        static bool compareFeh(const FehCurve &a, const double b)
        {
            return a.feh < b;
        }

        bool operator<(const struct FehCurve &b) const
        {
            return feh < b.feh;
        }

        double feh;
        std::vector<HeliumCurve> heliumCurves;
    };

  public:
    virtual ~MsRgbModel() {;}

/****************************************************************************************
 * Derive AGBt mass (actually the ZAMS mass for the appropriate AGBt star) for a given
 * metallicity, age, and helium abundance (if necessary).  Uses an array of pointers to functions.
 * Functions must return a double and have three double arguments (newFeH, newY, newLogAge)

 * Array indices are defined in evolve.h
****************************************************************************************/
    virtual Isochrone* deriveIsochrone(double, double, double) const = 0;

/****************************************************************************************
Derive WD precursor age for a given metallicity, calling in turn wd_prec_g_lage to
interpolate among Girardi or wd_prec_c_lage to interpolate among Chaboyer isochrones in
mass and age.

Distributed most of the code to the respective subroutines, leaving only those to be
modified for different model sets.
****************************************************************************************/
    virtual double wdPrecLogAge(double, double, double) const = 0;
    virtual void restrictToFilters(const std::vector<std::string>&) = 0;

    virtual double getMinAge() const { return ageLimit.first; }
    virtual double getMaxAge() const { return ageLimit.second; }

    virtual void setArtificialMinAge(double limit) { ageLimit.first  = limit; }
    virtual void setArtificialMaxAge(double limit) { ageLimit.second = limit; }

    virtual double getMinFeh() const { return fehCurves.front().feh; }
    virtual double getMaxFeh() const { return fehCurves.back().feh; }

    std::vector<std::string> getAvailableFilters() const { return availableFilters; }

  protected:
    virtual std::string getFileName (std::string) const = 0;

    std::pair<double, double> ageLimit;
    std::vector<std::string> availableFilters;

    std::vector<struct FehCurve> fehCurves;
};
#endif
