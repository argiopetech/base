#ifndef IO_STAR_PARAMS_MCMC_HPP
#define IO_STAR_PARAMS_MCMC_HPP

#include "BackingStore.hpp"
#include "Records.hpp"


class StarParams_FileBackingStore : public FileBackingStore<StarParamsRecord, double>
{
  public:
    StarParams_FileBackingStore(std::string);

    ~StarParams_FileBackingStore() override = default;

    void save(StarParamsRecord) override;

  private:
    void header(double) override { ; }
};

#endif
