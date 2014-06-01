#ifndef ALTHAUSWMODEL_HPP
#define ALTHAUSWMODEL_HPP

#include <string>
#include <utility>
#include <vector>

#include "CarbonlessWdModel.hpp"

class AlthausWdModel : public CarbonlessWdModel
{
  public:
    virtual void loadModel (std::string path);
};
#endif
