#ifndef GDSEDMAG_H
#define GDSEDMAG_H

#include <string>
#include <vector>

#include "GenericMsModel.hpp"

class DsedMsModel : public GenericMsModel
{
  public:
    DsedMsModel() {;}
    virtual ~DsedMsModel() {;}

    virtual bool isSupported(FilterSetName filterSet) const
        { return  filterSet == FilterSetName::UBVRIJHK || filterSet == FilterSetName::SDSS; }

  protected:
    virtual std::string getFileName (std::string path) const
        { return path + "dsed/dsed_old.model"; }
};

#endif
