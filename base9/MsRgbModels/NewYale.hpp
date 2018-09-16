#ifndef NEW_YALEMSMODEL_HPP
#define NEW_YALEMSMODEL_HPP

#include <string>
#include <vector>

#include "GenericMsModel.hpp"

class NewYaleMsModel : public GenericMsModel
{
  public:
    NewYaleMsModel() {;}
    virtual ~NewYaleMsModel() {;}

  protected:
    virtual std::string getFileName (std::string path) const
        { return path + "yale_2018/yale.model"; }
};

#endif
