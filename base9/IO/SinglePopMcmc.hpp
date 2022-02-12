#ifndef IO_SINGLE_POP_MCMC_HPP
#define IO_SINGLE_POP_MCMC_HPP

#include "BackingStore.hpp"
#include "Records.hpp"


class SinglePopMcmc_FileBackingStore : public FileBackingStore<SinglePopMcmcRecord, std::array<double, NPARAMS> const&>
{
  public:
    SinglePopMcmc_FileBackingStore(std::string);

    ~SinglePopMcmc_FileBackingStore() override = default;

    void save(SinglePopMcmcRecord) override;

  protected:
    void header(std::array<double, NPARAMS> const &params) override { header(params, false); }
    void header(std::array<double, NPARAMS> const&, bool);
};

#endif
