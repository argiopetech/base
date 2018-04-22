#ifndef IO_MULTI_POP_MCMC_HPP
#define IO_MULTI_POP_MCMC_HPP

#include "BackingStore.hpp"
#include "Records.hpp"


class MultiPopMcmc_FileBackingStore : public FileBackingStore<MultiPopMcmcRecord, std::array<double, NPARAMS> const&>
{
  public:
    MultiPopMcmc_FileBackingStore(std::string);
    
    ~MultiPopMcmc_FileBackingStore() = default;
    
    void save(Iteration, MultiPopMcmcRecord);

  protected:
    void header(std::array<double, NPARAMS> const&);
};

#endif
