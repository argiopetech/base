#include "McmcApplication.hpp"
#include "Model.hpp"
#include "Settings.hpp"

class MpiMcmcApplication : public McmcApplication
{
  public:
    MpiMcmcApplication(Settings &s)
        : settings(s), evoModels(makeModel(settings))
    {}

    virtual ~MpiMcmcApplication() {}

    virtual int run();

    void propClustBigSteps (Cluster &, const struct ifmrMcmcControl &) const;
    void propClustIndep (Cluster &, const struct ifmrMcmcControl &) const;
    void propClustCorrelated (Cluster &, const struct ifmrMcmcControl &) const;

  private:
    const Settings settings;
    const Model evoModels;
};
