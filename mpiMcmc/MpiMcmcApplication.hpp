#include "McmcApplication.hpp"

class MpiMcmcApplication : public McmcApplication
{
  public:
    virtual ~MpiMcmcApplication() {}

    virtual int run(int argc, char *argv[]);
};
