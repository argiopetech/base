#include <functional>
#include <random>

#include "mpiMcmc.hpp"
#include "McmcApplication.hpp"
#include "Model.hpp"
#include "Settings.hpp"
#include "Utility.hpp"

using std::hash;

class MpiMcmcApplication : public McmcApplication
{
  public:
    MpiMcmcApplication(Settings &s);

    virtual ~MpiMcmcApplication() {}

    virtual int run();

    Cluster propClustBigSteps (const Cluster&, const struct ifmrMcmcControl &, const std::array<double, NPARAMS> &);
    Cluster propClustIndep (Cluster, const struct ifmrMcmcControl &, const std::array<double, NPARAMS> &, double scale = 1.0);
    Cluster propClustCorrelated (Cluster, const struct ifmrMcmcControl &);

    double logPostStep (Cluster &, double, const std::vector<int>&);

  private:
    void scaleStepSizes (std::array<double, NPARAMS> &);

    double wdGridMass (int) const;

    const Model evoModels;
    const Settings settings;

    Cluster clust;
    vector<StellarSystem> systems;

    base::utility::ThreadPool pool;

    struct ifmrMcmcControl ctrl;

    int nSave = 100;

    const int trialIter;
};
