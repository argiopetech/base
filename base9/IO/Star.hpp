#ifndef IO_STAR_MCMC_HPP
#define IO_STAR_MCMC_HPP

#include "BackingStore.hpp"
#include "Records.hpp"


class Star_SqlBackingStore : public SqlBackingStore<StarRecord>
{
  public:
    Star_SqlBackingStore(const RunData&);
    Star_SqlBackingStore(const Star_SqlBackingStore&) = delete;
    Star_SqlBackingStore(std::string);

    ~Star_SqlBackingStore() override;

    void save(StarRecord) override;

  private:
    sqlite3_stmt *insertStar = nullptr;
    sqlite3_stmt *insertPhotometry = nullptr;

    void buildInsertStatement() override;
};

#endif
