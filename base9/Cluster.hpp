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
        stepSize.fill(0.0);
        mean.fill(0.0);
        parameter.fill(0.0);

        initModels (&evoModels);
    }
    ~Cluster() {;}

    void readClust (FILE * pFile);
    void writeClust (FILE * pFile);

    struct model evoModels;

    std::array<double, 3> betaAgeMod;
    std::array<double, NPARAMS> parameter;
    std::array<double, NPARAMS> stepSize;
    std::array<double, NPARAMS> mean;
    std::array<double, NPARAMS> priorVar;
    std::array<double, NPARAMS> priorMean;

    int nStars = 0;
    int photometrySet;
    double M_wd_up = 8.0;
    double betamabs = 0.0;
    double betaFabs = 0.0;
    double betaFY = 0.0;
    double AGBt_zmass = 0.0;
    double varScale = 1.0;
    double carbonicity = 0.0;
};

#endif

