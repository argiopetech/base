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

  private:
    const Settings settings;
    const Model evoModels;
};
