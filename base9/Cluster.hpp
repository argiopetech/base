#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <array>

#include "constants.hpp"
#include "Model.hpp"

class Cluster
{
  public:
    Cluster()
    {
        setM_wd_up(8.0);
    }

    void setParam(int, double);
    double getParam(int) const;

    void setM_wd_up(double);
    double getM_wd_up() const { return M_wd_up; }

    double getLogMassNorm() const { return logMassNorm; }

    double logPrior(const Model&) const;
    double logPriorMass(double) const;

    std::array<double, NPARAMS> priorVar;
    std::array<double, NPARAMS> priorMean;
    std::array<double, NPARAMS> mean;

    int photometrySet;
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

  private:
    double logMassNorm;
    double M_wd_up;
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
