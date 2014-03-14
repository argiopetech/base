#ifndef FILTERSET_HPP
#define FILTERSET_HPP

#include <array>
#include <string>

#include "constants.hpp"

class FilterSet
{
  public:
    FilterSet(std::array<std::string, FILTS> filts)
        : filterNames(filts)
    {}
    virtual ~FilterSet() {};

    std::string getFilterName (int index) const { return filterNames.at(index); };
    virtual std::array<double, FILTS> calcAbsCoeffs() const = 0;

  protected:
    const std::array<std::string, FILTS> filterNames;
};


class UBVRIJHK : public FilterSet
{
  public:
    UBVRIJHK();

    virtual std::array<double, FILTS> calcAbsCoeffs() const;
};

class SDSS : public FilterSet
{
  public:
    SDSS();

    virtual std::array<double, FILTS> calcAbsCoeffs() const;
};

class ACS : public FilterSet
{
  public:
    ACS();

    virtual std::array<double, FILTS> calcAbsCoeffs() const;
};

class UVIS : public FilterSet
{
  public:
    UVIS();

    virtual std::array<double, FILTS> calcAbsCoeffs() const;
};

#endif
