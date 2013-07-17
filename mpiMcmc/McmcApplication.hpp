#include <random>

class McmcApplication
{
  public:
    McmcApplication(int seed)
        : gen(seed)
    {}

    virtual ~McmcApplication() {}

    virtual int run() = 0;

  protected:
    // Decides whether to accept a proposed cluster property
    bool acceptClustMarg (const double logPostCurr, const double logPostProp);

    double acceptanceRatio() const { return double(accepted) / (accepted + rejected); }

    std::mt19937 gen;

  private:
    int accepted = 0;
    int rejected = 0;
};
