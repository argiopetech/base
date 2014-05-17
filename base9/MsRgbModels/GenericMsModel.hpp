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

    virtual void loadModel(std::string, FilterSetName);

    virtual double deriveAgbTipMass(double, double, double);
    virtual Isochrone deriveIsochrone(double, double, double) const;

    virtual std::vector<double> msRgbEvol (double) const;

    virtual double wdPrecLogAge(double, double, double) const;
    virtual void restrictToFilters(const std::vector<std::string>&);

  private:
    std::vector<std::string> availableFilters;

    Isochrone deriveIsochrone_oneY(const std::vector<int>&, double, double) const;
    Isochrone deriveIsochrone_manyY(const std::vector<int>&, double, double, double) const;

    double wdPrecLogAge_oneY(double, double) const;
    double wdPrecLogAge_manyY(double, double, double) const;
};
#endif
