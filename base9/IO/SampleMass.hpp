#ifndef IO_SAMPLE_MASS_HPP
#define IO_SAMPLE_MASS_HPP

#include "BackingStore.hpp"
#include "Records.hpp"


class SampleMass_FileBackingStore : public FileBackingStore<std::vector<SampleMassRecord>, std::vector<SampleMassRecord>>
{
  public:
    SampleMass_FileBackingStore(std::string);

    ~SampleMass_FileBackingStore() override = default;

    void save(std::vector<SampleMassRecord>) override;

  protected:
    void header(std::vector<SampleMassRecord>) override;

  private:
    // Default value to 11 to match the base::utility::format
    size_t longestStarIdLength = 11;
};

#endif
