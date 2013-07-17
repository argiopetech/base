#include <random>

#include "McmcApplication.hpp"
#include "Model.hpp"
#include "Settings.hpp"

class MpiMcmcApplication : public McmcApplication
{
  public:
    MpiMcmcApplication(Settings &s)
        : McmcApplication(s.seed), settings(s), evoModels(makeModel(settings))
    {}

    virtual ~MpiMcmcApplication() {}

    virtual int run();

    void propClustBigSteps (Cluster &, const struct ifmrMcmcControl &);
    void propClustIndep (Cluster &, const struct ifmrMcmcControl &);
    void propClustCorrelated (Cluster &, const struct ifmrMcmcControl &);

  private:
    const Settings settings;
    const Model evoModels;
};
