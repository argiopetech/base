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

    double logPostStep(Chain &, const std::array<double, N_WD_MASS1> &, Cluster &, double, const std::vector<int>&, std::array<double, FILTS>&, std::array<double, FILTS>&);

  private:
    const Model evoModels;
    const Settings settings;

    base::utility::ThreadPool pool;
};
