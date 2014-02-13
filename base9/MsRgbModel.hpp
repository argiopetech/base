#ifndef MSRGBEVOL_HPP
#define MSRGBEVOL_HPP

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "constants.hpp"
#include "structures.hpp"

class MsRgbModel
{
  public:
    virtual ~MsRgbModel() {;}

/****************************************************************************************
 * Derive AGBt mass (actually the ZAMS mass for the appropriate AGBt star) for a given
 * metallicity, age, and helium abundance (if necessary).  Uses an array of pointers to functions.
 * Functions must return a double and have three double arguments (newFeH, newY, newLogAge)

 * Array indices are defined in evolve.h
****************************************************************************************/
    virtual double deriveAgbTipMass(const std::vector<int>&, double, double, double) = 0;

/****************************************************************************************
Perform interpolation via  calls to getGirardiMags() or similar.
The former does 3-D interpolation of the Girardi isochrones.

deriveAgbTipMass() needs to be called first
****************************************************************************************/
    std::array<double, FILTS> msRgbEvol (const std::vector<int>&, double);


/****************************************************************************************
Derive WD precursor age for a given metallicity, calling in turn wd_prec_g_lage to
interpolate among Girardi or wd_prec_c_lage to interpolate among Chaboyer isochrones in
mass and age.

Distributed most of the code to the respective subroutines, leaving only those to be
modified for different model sets.
****************************************************************************************/
    virtual double wdPrecLogAge(double, double) = 0;

    virtual void loadModel(std::string, MsFilter) = 0;

    double getMinAge() const { return ageLimit.first; }
    double getMaxAge() const { return ageLimit.second; }

    const struct globalIso& getIsochrone() const { return isochrone; }

  protected:
    virtual int numFilts() const = 0;

    std::pair<double, double> ageLimit;

    struct globalIso isochrone;
};
#endif
