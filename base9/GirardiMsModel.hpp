#ifndef GGIRMAG_H
#define GGIRMAG_H

#include <string>

const int N_GIR_Z         = 8;    /* number of metallicities in Girardi isochrones */
const int N_GIR_AGES      = 50;   /* number of ages in Girardi isochonres */
const int N_GIR_FILTS     = 8;
const int MAX_GIR_ENTRIES = 10000;        /* max Girardi entries for given Z */

class GirardiMsModel : public MsRgbModel
{
  public:
    GirardiMsModel() {;}
    virtual ~GirardiMsModel() {;}

    virtual double deriveAgbTipMass(double, double, double);
    virtual double msRgbEvol(double);
    virtual double wdPrecLogAge(double, double, double);
    virtual void loadModel(std::string, MsFilterSet);

  private:
    double interpInMass (int, double, int, double *);
};

#endif
