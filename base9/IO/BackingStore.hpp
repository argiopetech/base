#ifndef IO_BACKING_STORE_HPP
#define IO_BACKING_STORE_HPP

#include <array>
#include <fstream>
#include <string>

// Facilitates other than busy looping on busy SQL databases
#include <chrono>
#include <thread>

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


// TODO: Make this follow NVI
// http://www.gotw.ca/publications/mill18.htm
template <typename T>
class BackingStore
{
  public:
    virtual ~BackingStore() = default;

    virtual void save(T) = 0;

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
    Iteration nextIteration() override;

  protected:
    std::shared_ptr<sqlite3> db;
    sqlite3_int64 run = -1;

    void execOnly(const std::string, const std::string);
    int execOnlyRet(const std::string);
    void execOnlyInTransaction(const std::string, const std::string);

    void dbErrorIf(int, const std::string);

    void dbStepAndReset(sqlite3_stmt*, const std::string);
    void dbPrepare(const std::string, sqlite3_stmt**, const std::string);

    void beginTransaction(const std::string);
    void endTransaction(const std::string);

  private:
    sqlite3_stmt *insertIteration = nullptr;

    void openDb(const std::string);
    void ensureTables();
    void generateRun();
    virtual void buildInsertStatement();

    // Other utilities
    void enableForeignKeys();
};

#include "SqlBackingStore.ipp"

#endif
