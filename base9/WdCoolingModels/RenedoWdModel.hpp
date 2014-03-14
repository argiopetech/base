#ifndef RENEDOWDMODEL_HPP
#define RENEDOWDMODEL_HPP

#include <string>
#include <utility>
#include <vector>

#include "CarbonlessWdModel.hpp"

class RenedoWdModel : public CarbonlessWdModel
{
  public:
    virtual void loadModel (std::string path, FilterSetName);
};
#endif
