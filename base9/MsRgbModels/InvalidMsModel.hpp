#ifndef INVALIDMSMODEL_HPP
#define INVALIDMSMODEL_HPP

#include <string>
#include <vector>

#include "../MsRgbModel.hpp"

class InvalidMsModel : public MsRgbModel, public InvalidModel
{
  public:
    InvalidMsModel() {;}
    virtual ~InvalidMsModel() {;}

    virtual void loadModel(std::string, FilterSetName) {;}
    virtual void restrictToFilters(const std::vector<std::string>&) {;}

    virtual double deriveAgbTipMass(double, double, double)
    {
        throw InvalidModelError("Called deriveAgbTipMass() in invalid MainSequence model");
    }

    virtual double wdPrecLogAge(double, double, double) const
    {
        throw InvalidModelError("Called wdPrecLogAge() in invalid MainSequence model");
    }

    virtual Isochrone deriveIsochrone(double, double, double) const
    {
        throw InvalidModelError("Called deriveIsochrone() in invalid MainSequence model");
    }

    virtual std::vector<double> msRgbEvol(double) const
    {
        throw InvalidModelError("Called msRgbEvol() in invalid MainSequence model");
    }

    virtual double getMinAge() const
    {
        throw InvalidModelError("Called getMinAge() in invalid MainSequence model");
    }

    virtual double getMaxAge() const
    {
        throw InvalidModelError("Called getMaxAge() in invalid MainSequence model");
    }

  protected:
    virtual std::string getFileName (std::string) const
    {
        throw InvalidModelError("Called getFileName() in invalid MainSequence model");
    }
};

#endif
