#ifndef GGIRMAG_H
#define GGIRMAG_H

#include <string>
#include <vector>

#include "GenericMsModel.hpp"

class GirardiMsModel : public GenericMsModel
{
  public:
    GirardiMsModel() {;}
    virtual ~GirardiMsModel() {;}

    virtual bool isSupported(FilterSetName filterSet) const
        { return filterSet == FilterSetName::UBVRIJHK || filterSet == FilterSetName::ACS; }

  protected:
    virtual std::string getFileName (std::string path) const
        { return path + "girardi/girardi.model"; }
};

#endif
