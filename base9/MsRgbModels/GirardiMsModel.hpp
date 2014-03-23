#ifndef GGIRMAG_H
#define GGIRMAG_H

#include <string>
#include <vector>

#include "../MsRgbModel.hpp"

const int N_GIR_Z         = 8;    /* number of metallicities in Girardi isochrones */
const int N_GIR_AGES      = 50;   /* number of ages in Girardi isochonres */
const int N_GIR_FILTS     = 8;
const int MAX_GIR_ENTRIES = 10000;        /* max Girardi entries for given Z */

class GirardiMsModel : public MsRgbModel
{
  public:
    GirardiMsModel() {;}
    virtual ~GirardiMsModel() {;}

    virtual void loadModel(std::string, FilterSetName);

    virtual double deriveAgbTipMass(const std::vector<int>&, double, double, double);
    virtual double wdPrecLogAge(double, double);

    virtual Isochrone deriveIsochrone(const std::vector<int>&, double, double, double) const;

    virtual bool isSupported(FilterSetName);

  protected:
    virtual int numFilts() const { return N_GIR_FILTS; }

  private:
    double modelZSolar = 0.019;
};

#endif
