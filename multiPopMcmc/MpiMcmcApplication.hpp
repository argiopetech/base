#ifndef MPIMCMCAPPLICATION_HPP
#define MPIMCMCAPPLICATION_HPP

#include <functional>
#include <random>
#include <tuple>

#include "IO/Records.hpp"
#include "mpiMcmc.hpp"
#include "McmcApplication.hpp"
#include "Model.hpp"
#include "Settings.hpp"
#include "Utility.hpp"


struct DualPopCluster
{
    std::array<Cluster, 2> clust;

    double lambda = 0.5;
};

class MpiMcmcApplication
{
  public:
    MpiMcmcApplication(Settings &s, MultiPopBackingStore*, StarParamsBackingStore*);

    virtual ~MpiMcmcApplication() {}

    virtual int run();

    DualPopCluster propClustBigSteps (DualPopCluster, const struct ifmrMcmcControl &, const std::array<double, NPARAMS> &);
    DualPopCluster propClustIndep (DualPopCluster, const struct ifmrMcmcControl &, const std::array<double, NPARAMS> &, double scale);
    DualPopCluster propClustCorrelated (DualPopCluster, const struct ifmrMcmcControl&, const Matrix<double, NPARAMS, NPARAMS>&);

    std::tuple<double, std::vector<double>> logPostStep (DualPopCluster &, double);

    void mainRun(std::function<DualPopCluster(DualPopCluster)> propose, std::function<double(DualPopCluster&)> logPost, std::ofstream &fout, std::array<double, NPARAMS> priorVar, DualPopCluster clust, int iters, int thin);


  private:
    void scaleStepSizes (std::array<double, NPARAMS> &, double);
    void allocateSSEMem();

    double wdGridMass (int) const;

    Model evoModels;
    const Settings settings;

    std::mt19937 gen;

    DualPopCluster clust;

    std::unique_ptr<MultiPopBackingStore> mcmcStore;
    std::unique_ptr<StarParamsBackingStore> paramsStore;

    std::vector<StellarSystem> systems;
    std::vector<StellarSystem> mainRunSystems;

    double* sysVars = nullptr;
    double* sysVar2 = nullptr;
    double* sysObs  = nullptr;

    size_t howManyFilts = 0;
    size_t howManyFiltsAligned = 0;
    base::utility::ThreadPool pool;

    struct ifmrMcmcControl ctrl;
};

#endif
