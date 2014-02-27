#ifndef FILTERSET_HPP
#define FILTERSET_HPP

#include <array>
#include <string>

#include "constants.hpp"

class MsFilterSet
{
  public:
    MsFilterSet(std::array<std::string, FILTS> filts)
        : filterNames(filts)
    {}
    virtual ~MsFilterSet() {};

    std::string getFilterName (int index) const { return filterNames.at(index); };
    virtual std::array<double, FILTS> calcAbsCoeffs() const = 0;

  protected:
    const std::array<std::string, FILTS> filterNames;
};


class UBVRIJHK : public MsFilterSet
{
  public:
    UBVRIJHK();

    virtual std::array<double, FILTS> calcAbsCoeffs() const;
};

class SDSS : public MsFilterSet
{
  public:
    SDSS();

    virtual std::array<double, FILTS> calcAbsCoeffs() const;
};

class ACS : public MsFilterSet
{
  public:
    ACS();

    virtual std::array<double, FILTS> calcAbsCoeffs() const;
};

class UVIS : public MsFilterSet
{
  public:
    UVIS();

    virtual std::array<double, FILTS> calcAbsCoeffs() const;
};

#endif
