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


class StarParams_SqlBackingStore : public SqlBackingStore<StarParamsRecord>
{
  public:
    StarParams_SqlBackingStore(const RunData&);
    StarParams_SqlBackingStore(const StarParams_SqlBackingStore&) = delete;
    StarParams_SqlBackingStore(std::string);

    ~StarParams_SqlBackingStore();

    void save(Iteration, StarParamsRecord);

  private:
    sqlite3_stmt *insertFsLike   = nullptr;
    sqlite3_stmt *insertStarData = nullptr;

    void buildInsertStatement();
};

#endif
