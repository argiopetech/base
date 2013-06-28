#ifndef CONSTANTS_H
#define CONSTANTS_H

const double EPSILON = 1.e-15;
const double EPS     = 0.000001;
const int FILTS      = 14;   // UBVRIJHK, respectively
const int MAG_LIST   = 26;   // list length for mag vs. S/N
const int MODEL_LIST = 118;  // age/Mv/V-I entries in WD model array

const int MSRG       = 1;    // Main-sequence or red giant
const int WD         = 3;    // White dwarf
const int NSBH       = 4;    // Neutron star or black hole
const int BD         = 5;    // Brown dwarf/planet
const int DNE        = 9;    // Does not exist (zams mass > 100 Msun or < 0.0001)

const double R_sun   = 6.96342e10;

const int NPARAMS            = 9;

// Filter sets
const int UBVRIJHK   = 0;
const int ACS        = 1;
const int SDSS       = 2;

// Define model sets
//Main sequence evolution
const int GIR        = 0;
const int CHABHELIUM = 1;
const int YALE       = 2;
const int DSED       = 3;

//IFMR
const int WEIDEMANN  = 0;
const int WILLIAMS   = 1;
const int SALARISLIN = 2;
const int SALARISPW  = 3;
const int LINEAR     = 4;

//WD Cooling
const int WOOD       = 0;
const int MONTGOMERY = 1;
const int ALTHAUS    = 2;
const int RENEDO     = 3;

//WD Atmosphere
const int BERGERON   = 0;

//Brown Dwarf/Planet
const int NONE       = 0;
const int BARAFFE    = 1;

const int ZAMS       = 0;
const int NOW        = 1;

const int LOW        = 0;
const int HIGH       = 1;

const int DA         = 0;
const int DB         = 1;

#endif
