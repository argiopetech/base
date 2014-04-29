#ifndef MSRGBEVOL_HPP
#define MSRGBEVOL_HPP

#include <array>
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

    virtual void loadModel(std::string, FilterSetName);

/****************************************************************************************
 * Derive AGBt mass (actually the ZAMS mass for the appropriate AGBt star) for a given
 * metallicity, age, and helium abundance (if necessary).  Uses an array of pointers to functions.
 * Functions must return a double and have three double arguments (newFeH, newY, newLogAge)

 * Array indices are defined in evolve.h
****************************************************************************************/
    virtual double deriveAgbTipMass(const std::vector<int>&, double, double, double);
    virtual Isochrone deriveIsochrone(const std::vector<int>&, double, double, double) const;

/****************************************************************************************
Perform interpolation via  calls to getGirardiMags() or similar.
The former does 3-D interpolation of the Girardi isochrones.

deriveAgbTipMass() needs to be called first
****************************************************************************************/
    std::array<double, FILTS> msRgbEvol (const std::vector<int>&, double) const;

/****************************************************************************************
Derive WD precursor age for a given metallicity, calling in turn wd_prec_g_lage to
interpolate among Girardi or wd_prec_c_lage to interpolate among Chaboyer isochrones in
mass and age.

Distributed most of the code to the respective subroutines, leaving only those to be
modified for different model sets.
****************************************************************************************/
    virtual double wdPrecLogAge(double, double);

    double getMinAge() const { return ageLimit.first; }
    double getMaxAge() const { return ageLimit.second; }

    const Isochrone& getIsochrone() const { return isochrone; }

  protected:
    virtual int numFilts() const = 0;
    virtual std::string getFileName (std::string) const = 0;

    std::pair<double, double> ageLimit;

    Isochrone isochrone;
    std::vector<struct FehCurve> fehCurves;

  private:
    std::vector<std::string> availableFilters;
};
#endif
