#ifndef MONTGOMERYWDMODEL_HPP
#define MONTGOMERYWDMODEL_HPP

#include <string>
#include <utility>
#include <vector>

#include "../WdCoolingModel.hpp"

class MontgomeryWdModel : public WdCoolingModel
{
  public:
    virtual void loadModel (std::string path, FilterSetName);
    virtual std::pair<double, double> wdMassToTeffAndRadius (double logAge, double x_carbon, double wdPrecLogAge, double wdMass) const;
};
#endif
