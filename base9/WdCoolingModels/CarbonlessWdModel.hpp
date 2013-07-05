#ifndef CARBONLESSWMODEL_HPP
#define CARBONLESSWMODEL_HPP

#include <string>
#include <utility>
#include <vector>

#include "../WdCoolingModel.hpp"

class CarbonlessWdModel : public WdCoolingModel
{
  public:
    virtual std::pair<double, double> wdMassToTeffAndRadius (double logAge, double x_carbon, double wdPrecLogAge, double wdMass) const;
};
#endif
