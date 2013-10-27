#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <array>

#include "constants.hpp"

class Cluster
{
  public:
    double getParam(int) const;
    void setParam(int, double);

    std::array<double, 3> betaAgeMod;
    std::array<double, NPARAMS> priorVar;
    std::array<double, NPARAMS> priorMean;
    std::array<double, NPARAMS> mean;

    int photometrySet;
    double M_wd_up = 8.0;
    double AGBt_zmass = 0.0;
    double varScale = 1.0;

    double age = 0.0;
    double yyy = 0.0;
    double feh = 0.0;
    double mod = 0.0;
    double abs = 0.0;
    double carbonicity = 0.0;
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
