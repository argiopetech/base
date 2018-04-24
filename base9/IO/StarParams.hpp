#ifndef IO_STAR_PARAMS_MCMC_HPP
#define IO_STAR_PARAMS_MCMC_HPP

#include "BackingStore.hpp"
#include "Records.hpp"


class StarParams_FileBackingStore : public FileBackingStore<StarParamsRecord, double>
{
  public:
    StarParams_FileBackingStore(std::string);
    
    ~StarParams_FileBackingStore() = default;
    
    void save(Iteration, StarParamsRecord);

  private:
    void header(double);
};

#endif
