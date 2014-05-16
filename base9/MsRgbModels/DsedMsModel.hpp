#ifndef GDSEDMAG_H
#define GDSEDMAG_H

#include <string>
#include <vector>

#include "../MsRgbModel.hpp"

const int N_DSED_FILTS = 13;

class DsedMsModel : public MsRgbModel
{
  public:
    DsedMsModel() {;}
    virtual ~DsedMsModel() {;}

    virtual bool isSupported(FilterSetName filterSet) const
        { return  filterSet == FilterSetName::UBVRIJHK || filterSet == FilterSetName::SDSS; }

  protected:
    virtual int numFilts() const
        { return N_DSED_FILTS; }

    virtual std::string getFileName (std::string path) const
        { return path + "dsed/dsed_old.model"; }
};

#endif
