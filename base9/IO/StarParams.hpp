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


class StarParams_SqlBackingStore : public SqlBackingStore<StarParamsRecord>
{
  public:
    StarParams_SqlBackingStore(const RunData&);
    StarParams_SqlBackingStore(const SqlBackingStore&);
    StarParams_SqlBackingStore(const StarParams_SqlBackingStore&) = delete;
    StarParams_SqlBackingStore(std::string);

    ~StarParams_SqlBackingStore() override;

    void save(StarParamsRecord) override;

  private:
    sqlite3_stmt *insert = nullptr;

    void buildInsertStatement() override;
};

#endif
