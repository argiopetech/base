#ifndef PARSECMSMODEL_HPP
#define PARSECMSMODEL_HPP

#include <string>
#include <vector>

#include "GenericMsModel.hpp"

class ParsecMsModel : public GenericMsModel
{
  public:
    ParsecMsModel() {;}
    virtual ~ParsecMsModel() {;}

  protected:
    virtual std::string getFileName (std::string path) const
        { return path + "PanSTARRS/PanSTARRS.model"; }
};

#endif
