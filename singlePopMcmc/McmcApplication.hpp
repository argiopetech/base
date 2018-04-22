#ifndef MCMCAPPLICATION_HPP
#define MCMCAPPLICATION_HPP

#include <iostream>
#include <random>


class McmcApplication
{
  public:
    McmcApplication(uint32_t seed)
        : gen(uint32_t(seed * uint32_t(2654435761))) // Applies Knuth's multiplicative hash for obfuscation (TAOCP Vol. 3)
    {
        const int warmupIter = 10000;

        // std::cout << "Warming up generator..." << std::flush;

        gen.discard(warmupIter);

        // std::cout << " Done." << std::endl;
        // std::cout << "Generated " << warmupIter << " values." << std::endl;
    }

    virtual ~McmcApplication() {}

    double acceptanceRatio() const { return double(accepted) / (accepted + rejected); }

    void resetRatio()
    {
        accepted = 0;
        rejected = 0;
    }

  protected:
    // Decides whether to accept a proposed cluster property
    bool acceptClustMarg (const double logPostCurr, const double logPostProp);

    std::mt19937 gen;

  private:
    int accepted = 0;
    int rejected = 0;
};

#endif
