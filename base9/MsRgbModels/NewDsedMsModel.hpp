#ifndef NEWDSEDMSMODEL_HPP
#define NEWDSEDMSMODEL_HPP

#include <string>
#include <vector>

#include "GenericMsModel.hpp"

class NewDsedMsModel : public GenericMsModel
{
  public:
    NewDsedMsModel() {;}
    virtual ~NewDsedMsModel() {;}

    virtual bool isSupported(FilterSetName filterSet) const
        { return  filterSet == FilterSetName::UBVRIJHK || filterSet == FilterSetName::ACS || filterSet == FilterSetName::UVIS; }

  protected:
    virtual std::string getFileName (std::string path) const
        { return path + "dsed/dsed_new.model"; }
};

#endif
