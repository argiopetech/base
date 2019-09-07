#ifndef INVALIDMSMODEL_HPP
#define INVALIDMSMODEL_HPP

#include <string>
#include <numeric>
#include <vector>

#include "../MsRgbModel.hpp"

class InvalidMsModel : public MsRgbModel, public InvalidModel
{
  public:
    InvalidMsModel() {;}
    virtual ~InvalidMsModel() {;}

    virtual void loadModel(std::string) {;}
    virtual void restrictToFilters(const std::vector<std::string>&, bool) {;}

    virtual double wdPrecLogAge(double, double, double) const
    {
        throw InvalidModelError("Called wdPrecLogAge() in invalid MainSequence model");
    }

    virtual Isochrone* deriveIsochrone(double, double, double) const
    {
        throw InvalidModelError("Called deriveIsochrone() in invalid MainSequence model");
    }

    virtual double getMinAge() const
    {
        return std::numeric_limits<double>::lowest();
    }

    virtual double getMaxAge() const
    {
        return std::numeric_limits<double>::max();
    }

  protected:
    virtual std::string getFileName (std::string) const
    {
        throw InvalidModelError("Called getFileName() in invalid MainSequence model");
    }
};

#endif
