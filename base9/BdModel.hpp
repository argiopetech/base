#ifndef BDMODEL_HPP
#define BDMODEL_HPP

#include <array>
#include <string>
#include <vector>

#include "BdModel.hpp"

class BdModel
{
  public:
    virtual ~BdModel;

    virtual void loadModel (std::string) = 0;
    virtual void getMags (const std::vector<int>&, std::array<double, FILTS>&, double, double) = 0;
};
#endif
