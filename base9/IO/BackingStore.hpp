#ifndef IO_BACKING_STORE_HPP
#define IO_BACKING_STORE_HPP

#include <array>
#include <fstream>
#include <string>

#include "Cluster.hpp"
#include "constants.hpp"

#include "sqlite/sqlite3.h"


enum class AdaptiveMcmcStage { FixedBurnin, AdaptiveBurnin, AdaptiveMainRun, MainRun };

struct RunData
{
    const std::shared_ptr<sqlite3> db;
    const sqlite3_int64 run;
};


struct Iteration
{
    const int val;
};


template <typename T>
class BackingStore
{
  public:
    virtual ~BackingStore() = default;

    virtual void save(Iteration, T) = 0;

    virtual Iteration nextIteration() { return {++iteration}; }

  private:
    int iteration = 0;
};


template <typename T, typename HeaderData>
class FileBackingStore : public BackingStore<T>
{
  public:
    FileBackingStore(const std::string);

    virtual ~FileBackingStore();

  protected:
    std::ofstream fout;

    virtual void header(HeaderData) = 0;
};

#include "FileBackingStore.ipp"


template <typename T>
class SqlBackingStore : public BackingStore<T>
{
  public:
    SqlBackingStore(const std::string);
    SqlBackingStore(const RunData&);
    SqlBackingStore(const SqlBackingStore&);

    virtual ~SqlBackingStore();

    RunData runData() { return {db, run}; }
    Iteration nextIteration();

  protected:
    std::shared_ptr<sqlite3> db;
    sqlite3_int64 run = -1;

    int execOnly(const std::string);

    void dbErrorIf(int, const std::string);
    void dbStepAndReset(sqlite3_stmt*, const std::string);
    void dbPrepare(const std::string, sqlite3_stmt**, const std::string);

  private:
    sqlite3_stmt *insertIteration = nullptr;

    void openDb(const std::string);
    void ensureTables();
    void generateRun();
    void buildInsertStatement();
};

#include "SqlBackingStore.ipp"

#endif
