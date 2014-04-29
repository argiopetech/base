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

    virtual double deriveAgbTipMass(const std::vector<int>&, double, double, double)
    {
        throw InvalidModelError("Called deriveAgbTipMass() in invalid MainSequence model");
    }

    virtual double wdPrecLogAge(double, double)
    {
        throw InvalidModelError("Called wdPrecLogAge() in invalid MainSequence model");
    }

    virtual Isochrone deriveIsochrone(const std::vector<int>&, double, double, double) const
    {
        throw InvalidModelError("Called deriveIsochrone() in invalid MainSequence model");
    }

  protected:
    virtual int numFilts() const
    {
        throw InvalidModelError("Called numFilts() in invalid MainSequence model");
    }

    virtual std::string getFileName (std::string) const
    {
        throw InvalidModelError("Called getFileName() in invalid MainSequence model");
    }
};

#endif
