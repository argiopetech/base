#ifndef MPIMCMCAPPLICATION_HPP
#define MPIMCMCAPPLICATION_HPP

#include <functional>
#include <random>

#include "mpiMcmc.hpp"
#include "McmcApplication.hpp"
#include "Model.hpp"
#include "Settings.hpp"
#include "Utility.hpp"

class MpiMcmcApplication
{
  public:
    MpiMcmcApplication(Settings &s);

    virtual ~MpiMcmcApplication() {}

    virtual int run();

    Cluster propClustBigSteps (Cluster, const struct ifmrMcmcControl &, const std::array<double, NPARAMS> &);
    Cluster propClustIndep (Cluster, const struct ifmrMcmcControl &, const std::array<double, NPARAMS> &, double scale);
    Cluster propClustCorrelated (Cluster, const struct ifmrMcmcControl&, const Matrix<double, NPARAMS, NPARAMS>&);

    double logPostStep (Cluster &, double);

    void mainRun(std::function<Cluster(Cluster)> propose, std::function<double(Cluster&)> logPost, std::ofstream &fout, std::array<double, NPARAMS> priorVar, Cluster clust, int iters, int thin);


  private:
    void scaleStepSizes (std::array<double, NPARAMS> &, double);

    double wdGridMass (int) const;

    Model evoModels;
    const Settings settings;

    std::mt19937 gen;

    Cluster clust;
    Cluster mainClust;
    std::vector<StellarSystem> msSystems;
    std::vector<StellarSystem> wdSystems;

    double* sysVars = nullptr;
    double* sysVar2 = nullptr;
    double* sysObs  = nullptr;

    size_t howManyFilts = 0;
    size_t howManyFiltsAligned = 0;

    base::utility::ThreadPool pool;

    struct ifmrMcmcControl ctrl;
};

#endif
