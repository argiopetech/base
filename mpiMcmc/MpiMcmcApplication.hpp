#include <random>

#include "mpiMcmc.hpp"
#include "McmcApplication.hpp"
#include "Model.hpp"
#include "Settings.hpp"
#include "Utility.hpp"

class MpiMcmcApplication : public McmcApplication
{
  public:
    MpiMcmcApplication(Settings &s)
        : McmcApplication(s.seed), evoModels(makeModel(s)), settings(s), pool(2)
    {}

    virtual ~MpiMcmcApplication() {}

    virtual int run();

    void propClustBigSteps (Cluster &, const struct ifmrMcmcControl &);
    void propClustIndep (Cluster &, const struct ifmrMcmcControl &);
    void propClustCorrelated (Cluster &, const struct ifmrMcmcControl &);

    double logPostStep(Chain &, Cluster &, double, const std::vector<int>&, std::array<double, FILTS>&, std::array<double, FILTS>&);

  private:
    std::pair<std::array<double, N_MS_MASS1 * N_MS_MASS_RATIO>, array<double, N_MS_MASS1 * N_MS_MASS_RATIO>> initMsMassRatioGrid (const Chain&);
    std::array<double, N_WD_MASS1> initWdMass1Grid (const Chain&);

    double wdGridMass (int) const;

    const Model evoModels;
    const Settings settings;

    base::utility::ThreadPool pool;
};
