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
enum class MsFilter {UBVRIJHK, ACS, SDSS};

// Define model sets
//Main sequence evolution
enum class MsModel {GIRARDI, CHABHELIUM, YALE, DSED};

//IFMR
const int WEIDEMANN  = 0;
const int WILLIAMS   = 1;
const int SALARISLIN = 2;
const int SALARISPW  = 3;
const int LINEAR     = 4;

//WD Cooling
enum class WdModel {WOOD, MONTGOMERY, ALTHAUS, RENEDO};

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

//Be careful adding new sample types.  Their order is important.
//There are a few places in the code that test for SAMPLE_TYPE > or <.
const int MASS              = -2;
const int AGE_DURING_WANDER = -1;  // age sampling during AGE_WANDER (defined in samplers.h)

const int AGE               =  0;  // age sampling
const int YYY               =  1;  // helium sampling
const int FEH               =  2;  // metallicity sampling
const int MOD               =  3;  // modulus sampling
const int ABS               =  4;  // absorption sampling;
const int IFMR_INTERCEPT    =  5;
const int IFMR_SLOPE        =  6;
const int IFMR_QUADCOEF     =  7;

const int BERG_N_DA_LOG_G = 6;
const int BERG_N_DA_TEFF  = 57;
const int BERG_N_DB_LOG_G = 5;
const int BERG_N_DB_TEFF  = 31;
const int BERG_NFILTS     = 8;

#endif
