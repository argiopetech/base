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

    Cluster propClustBigSteps (const Cluster&, const struct ifmrMcmcControl &);
    Cluster propClustIndep (Cluster, const struct ifmrMcmcControl &, double scale = 1.0);
    Cluster propClustCorrelated (Cluster, const struct ifmrMcmcControl &);

    double logPostStep(const std::vector<Star> &, Cluster &, double, const std::vector<int>&, std::array<double, FILTS>&, std::array<double, FILTS>&);

  private:
//    std::pair<std::array<double, N_MS_MASS1 * N_MS_MASS_RATIO>, array<double, N_MS_MASS1 * N_MS_MASS_RATIO>> initMsMassRatioGrid (const Chain&);
//    std::array<double, N_WD_MASS1> initWdMass1Grid (const Chain&);

    double wdGridMass (int) const;

    const Model evoModels;
    const Settings settings;

    Cluster clust;
    vector<Star> stars;

    base::utility::ThreadPool pool;

    struct ifmrMcmcControl ctrl;

    int nSave = 100;
};
