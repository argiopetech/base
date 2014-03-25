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

    Cluster propClustBigSteps (Cluster, const struct ifmrMcmcControl &, const std::array<double, NPARAMS> &);
    Cluster propClustIndep (Cluster, const struct ifmrMcmcControl &, const std::array<double, NPARAMS> &, double scale);
    Cluster propClustCorrelated (Cluster, const struct ifmrMcmcControl&, const Matrix<double, NPARAMS, NPARAMS>&);

    double logPostStep (Cluster &, double, const std::vector<int>&);

    void mainRun(std::function<Cluster(Cluster)> propose, std::function<double(Cluster&)> logPost, ofstream &fout, std::array<double, NPARAMS> priorVar, Cluster clust, int iters, int thin);


  private:
    void scaleStepSizes (std::array<double, NPARAMS> &, double);

    double wdGridMass (int) const;

    const Model evoModels;
    const Settings settings;

    Cluster clust;
    Cluster mainClust;
    vector<StellarSystem> systems;

    base::utility::ThreadPool pool;

    struct ifmrMcmcControl ctrl;
};
