#ifndef IO_MULTI_POP_MCMC_HPP
#define IO_MULTI_POP_MCMC_HPP

#include "BackingStore.hpp"
#include "Records.hpp"


class MultiPopMcmc_FileBackingStore : public FileBackingStore<MultiPopMcmcRecord, std::array<double, NPARAMS> const&>
{
  public:
    MultiPopMcmc_FileBackingStore(std::string);

    ~MultiPopMcmc_FileBackingStore() override = default;

    void save(MultiPopMcmcRecord) override;

  protected:
    void header(std::array<double, NPARAMS> const &params) override { header(params, false); }
    void header(std::array<double, NPARAMS> const&, bool);
};


class MultiPopMcmc_SqlBackingStore : public SqlBackingStore<MultiPopMcmcRecord>
{
  public:
    MultiPopMcmc_SqlBackingStore(const RunData&);
    MultiPopMcmc_SqlBackingStore(const MultiPopMcmc_SqlBackingStore&) = delete;
    MultiPopMcmc_SqlBackingStore(std::string);

    ~MultiPopMcmc_SqlBackingStore() override;

    void save(MultiPopMcmcRecord) override;

  private:
    sqlite3_stmt *insert = nullptr;

    void buildInsertStatement() override;
};

#endif
