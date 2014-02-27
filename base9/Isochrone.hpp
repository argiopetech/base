#ifndef ISOCHRONE_HPP
#define ISOCHRONE_HPP

#include <array>
#include <limits>
#include <vector>

struct EvolutionaryPoint
{
    EvolutionaryPoint(int eep, double mass, std::array<double, FILTS> mags)
        : eep(eep), mass(mass), mags(mags)
    {;}

    ~EvolutionaryPoint()
    {;}

    // These are for use with search/sort functions in e.g., <algorithms>
    static bool compareEep(const EvolutionaryPoint &a, const double b)
    {
        return a.eep < b;
    }

    static bool compareMass(const EvolutionaryPoint &a, const double b)
    {
        return a.mass < b;
    }

    bool operator<(const struct EvolutionaryPoint &b) const
    {
        return mass < b.mass;
    }    

    int eep;
    double mass;

    std::array<double, FILTS> mags;
};

struct Isochrone
{
    Isochrone()
    {
        logAge = std::numeric_limits<double>::lowest();
    }

    Isochrone(double logAge, std::vector<EvolutionaryPoint> eeps)
        : logAge(logAge), eeps(eeps)
    {;}

    ~Isochrone()
    {;}

    static bool compareAgbTip(const Isochrone &a, const double b)
    {
        return a.agbTipMass() < b;
    }

    bool operator<(const struct Isochrone &b) const
    {
        return logAge < b.logAge;
    }    

    double agbTipMass() const
    {
        return eeps.back().mass;
    }

    double logAge;

    std::vector<EvolutionaryPoint> eeps;
};

#endif
