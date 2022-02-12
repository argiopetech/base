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

#endif
