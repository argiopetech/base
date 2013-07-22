#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <array>

#include "constants.hpp"

class Cluster
{
  public:
    Cluster()
    {
//        mean.fill(0.0);

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

    double getParam(int) const;
    void setParam(int, double);

    std::array<double, 3> betaAgeMod;
    std::array<double, NPARAMS> stepSize;
    std::array<double, NPARAMS> priorVar;
    std::array<double, NPARAMS> priorMean;

    int photometrySet;
    double M_wd_up = 8.0;
    double AGBt_zmass = 0.0;
    double varScale = 1.0;
    double carbonicity = 0.38; // Good default value, per Mike Montgomery

    double age = 0.0;
    double yyy = 0.0;
    double feh = 0.0;
    double mod = 0.0;
    double abs = 0.0;
    double ifmrIntercept = 0.0;
    double ifmrSlope = 0.0;
    double ifmrQuadCoef = 0.0;

};

class InvalidCluster : public std::range_error
{
  public:
    explicit InvalidCluster (const std::string& what_arg)
        : std::range_error(what_arg) {}

    explicit InvalidCluster (const char* what_arg)
        : std::range_error(what_arg) {}
};
#endif
