#ifndef IO_FIELD_STAR_LIKELIHOOD_MCMC_HPP
#define IO_FIELD_STAR_LIKELIHOOD_MCMC_HPP

#include "BackingStore.hpp"
#include "Records.hpp"


class FieldStarLikelihood_FileBackingStore : public FileBackingStore<FieldStarLikelihoodRecord, double>
{
  public:
    FieldStarLikelihood_FileBackingStore(std::string);

    ~FieldStarLikelihood_FileBackingStore() override = default;

    void save(FieldStarLikelihoodRecord) override;

  private:
    void header(double) override;
};


class FieldStarLikelihood_SqlBackingStore : public SqlBackingStore<FieldStarLikelihoodRecord>
{
  public:
    FieldStarLikelihood_SqlBackingStore(const RunData&);
    FieldStarLikelihood_SqlBackingStore(const FieldStarLikelihood_SqlBackingStore&) = delete;
    FieldStarLikelihood_SqlBackingStore(std::string);

    ~FieldStarLikelihood_SqlBackingStore() override;

    void save(FieldStarLikelihoodRecord) override;

  private:
    sqlite3_stmt *insert = nullptr;

    void buildInsertStatement() override;
};

#endif
