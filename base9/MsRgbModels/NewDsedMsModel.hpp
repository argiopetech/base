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

  protected:
    virtual std::string getFileName (std::string path) const
        { return path + "dsed/dsed_new.model"; }
};

#endif
