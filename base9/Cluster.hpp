#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <array>

#include "constants.hpp"

class Cluster
{
  public:
    Cluster()
    {
        mean.fill(0.0);
        parameter.fill(0.0);

        stepSize[AGE] = 0.005;
        stepSize[FEH] = 0.005;
        stepSize[MOD] = 0.005;
        stepSize[ABS] = 0.002;
        stepSize[YYY] = 0.002;
        stepSize[IFMR_INTERCEPT] = 0.01;
        stepSize[IFMR_SLOPE] = 0.008;
        stepSize[IFMR_QUADCOEF] = 0.008;
    }
    ~Cluster() {;}

    double getAge() const;
    double getY() const;
    double getFeH() const;
    double getMod() const;
    double getAbs() const;

    void setAge(const double newAge);
    void setY(const double newY);
    void setFeH(const double newFeH);
    void setMod(const double newMod);
    void setAbs(const double newAbs);

    std::array<double, 3> betaAgeMod;
    std::array<double, NPARAMS> parameter;
    std::array<double, NPARAMS> stepSize;
    std::array<double, NPARAMS> mean;
    std::array<double, NPARAMS> priorVar;
    std::array<double, NPARAMS> priorMean;

    int photometrySet;
    double M_wd_up = 8.0;
    double AGBt_zmass = 0.0;
    double varScale = 1.0;
    double carbonicity = 0.38; // Good default value, per Mike Montgomery
};

#endif
