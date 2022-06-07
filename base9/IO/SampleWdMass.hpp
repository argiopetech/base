#ifndef IO_SAMPLE_WD_MASS_HPP
#define IO_SAMPLE_WD_MASS_HPP

#include "BackingStore.hpp"
#include "Records.hpp"


class SampleWdMass_FileBackingStore : public FileBackingStore<std::vector<SampleWdMassRecord>, std::vector<SampleWdMassRecord>>
{
  public:
    SampleWdMass_FileBackingStore(std::string);

    ~SampleWdMass_FileBackingStore() override = default;

    void save(std::vector<SampleWdMassRecord>) override;

  protected:
    void header(std::vector<SampleWdMassRecord>) override;

  private:
    // Default value to 11 to match the base::utility::format
    size_t longestStarIdLength = 11;
};

#endif
