#ifndef MSRGBEVOL_HPP
#define MSRGBEVOL_HPP

#include <string>

#include "constants.hpp"

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
    virtual double deriveAgbTipMass(double, double, double) = 0;

/****************************************************************************************
Perform interpolation via  calls to getGirardiMags() or similar.
The former does 3-D interpolation of the Girardi isochrones.

deriveAgbTipMass() needs to be called first
****************************************************************************************/
    virtual double msRgbEvol(double) = 0;

/****************************************************************************************
Derive WD precursor age for a given metallicity, calling in turn wd_prec_g_lage to
interpolate among Girardi or wd_prec_c_lage to interpolate among Chaboyer isochrones in
mass and age.

Distributed most of the code to the respective subroutines, leaving only those to be
modified for different model sets.
****************************************************************************************/
    virtual double wdPrecLogAge(double, double, double) = 0;

    virtual void loadModel(std::string, MsFilterSet) = 0;
};
#endif
