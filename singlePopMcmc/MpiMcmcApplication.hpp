#ifndef MPIMCMCAPPLICATION_HPP
#define MPIMCMCAPPLICATION_HPP

#include <functional>
#include <random>

#include "Chain.hpp"
#include "IO/BackingStore.hpp"
#include "IO/Records.hpp"
#include "mpiMcmc.hpp"
#include "McmcApplication.hpp"
#include "Model.hpp"
#include "Settings.hpp"
#include "Utility.hpp"


class MpiMcmcApplication
{
  public:
    MpiMcmcApplication(Settings &s, SinglePopBackingStore*);

    virtual ~MpiMcmcApplication() {}

    virtual int run();

    Cluster propClustBigSteps (Cluster, const struct ifmrMcmcControl &, const std::array<double, NPARAMS> &);
    Cluster propClustIndep (Cluster, const struct ifmrMcmcControl &, const std::array<double, NPARAMS> &, double scale);
    Cluster propClustCorrelated (Cluster, const struct ifmrMcmcControl&, const Matrix<double, NPARAMS, NPARAMS>&);

    double logPostStep (Cluster &, double);

  private:
    void scaleStepSizes (std::array<double, NPARAMS> &, double);
    void allocateSSEMem();

    double wdGridMass (int) const;

    Model evoModels;
    const Settings settings;

    std::mt19937 gen;

    Cluster clust;

    std::unique_ptr<SinglePopBackingStore> backingStore;

    std::vector<StellarSystem> msSystems;
    std::vector<StellarSystem> msMainRun;

    std::vector<StellarSystem> wdSystems;
    std::vector<StellarSystem> wdMainRun;

    double* sysVars = nullptr;
    double* sysVar2 = nullptr;
    double* sysObs  = nullptr;

    size_t howManyFilts = 0;
    size_t howManyFiltsAligned = 0;

    base::utility::ThreadPool pool;

    struct ifmrMcmcControl ctrl;
};

#endif
