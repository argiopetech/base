#ifndef NEW_MONTGOMERYWDMODEL_HPP
#define NEW_MONTGOMERYWDMODEL_HPP

#include <string>
#include <utility>
#include <vector>

#include "../WdCoolingModel.hpp"

class NewMontgomeryWdModel : public WdCoolingModel
{
  public:
    virtual void loadModel (std::string path);
    virtual std::pair<double, double> wdMassToTeffAndRadius (double logAge, double x_carbon, double wdPrecLogAge, double wdMass) const;
};
#endif
