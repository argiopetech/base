#ifndef WOODWDMODEL_HPP
#define WOODWDMODEL_HPP

#include <string>
#include <utility>
#include <vector>

#include "CarbonlessWdModel.hpp"

class WoodWdModel : public CarbonlessWdModel
{
  public:
    virtual void loadModel (std::string path);
};
#endif
