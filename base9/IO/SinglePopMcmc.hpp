#ifndef IO_SINGLE_POP_MCMC_HPP
#define IO_SINGLE_POP_MCMC_HPP

#include "BackingStore.hpp"
#include "Records.hpp"


class SinglePopMcmc_FileBackingStore : public FileBackingStore<SinglePopMcmcRecord, std::array<double, NPARAMS> const&>
{
  public:
    SinglePopMcmc_FileBackingStore(std::string);
    
    ~SinglePopMcmc_FileBackingStore() = default;
    
    void save(Iteration, SinglePopMcmcRecord);

  protected:
    void header(std::array<double, NPARAMS> const&);
};

#endif
