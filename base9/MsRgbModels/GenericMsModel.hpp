#ifndef GENERICMSMODEL_HPP
#define GENERICMSMODEL_HPP

#include <string>
#include <utility>
#include <vector>

#include "../MsRgbModel.hpp"

class GenericMsModel : public MsRgbModel
{
  public:
    virtual ~GenericMsModel() {;}

    virtual void loadModel(std::string);

    virtual Isochrone* deriveIsochrone(double, double, double) const;

    virtual double wdPrecLogAge(double, double, double) const;
    virtual void restrictToFilters(const std::vector<std::string>&);

  private:
    Isochrone* deriveIsochrone_oneY(double, double) const;
    Isochrone* deriveIsochrone_manyY(double, double, double) const;
};
#endif
