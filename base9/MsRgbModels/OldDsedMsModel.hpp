#ifndef OLDDSEDMSMODEL_HPP
#define OLDDSEDMSMODEL_HPP

#include <string>
#include <vector>

#include "GenericMsModel.hpp"

class OldDsedMsModel : public GenericMsModel
{
  public:
    OldDsedMsModel() {;}
    virtual ~OldDsedMsModel() {;}

    virtual bool isSupported(FilterSetName filterSet) const
        { return  filterSet == FilterSetName::UBVRIJHK || filterSet == FilterSetName::SDSS; }

  protected:
    virtual std::string getFileName (std::string path) const
        { return path + "dsed/dsed_old.model"; }
};

#endif
