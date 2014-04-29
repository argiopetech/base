#ifndef GDSEDMAG_H
#define GDSEDMAG_H

#include <string>
#include <vector>

#include "../MsRgbModel.hpp"

const int N_DSED_Z         = 9;    /* number of metallicities in Chaboyer-Dotter isochrones */
const int N_DSED_AGES      = 52;   /* number of ages in Chaboyer-Dotter isochonres */
const int N_DSED_FILTS     = 8;
const int MAX_DSED_ENTRIES = 370;

class DsedMsModel : public MsRgbModel
{
  public:
    DsedMsModel() {;}
    virtual ~DsedMsModel() {;}

    virtual bool isSupported(FilterSetName) const;

  protected:
    virtual int numFilts() const { return N_DSED_FILTS; }
    virtual std::string getFileName (std::string) const;
};

#endif
