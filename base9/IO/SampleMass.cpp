#include <array>
#include <string>

#include "SampleMass.hpp"
#include "Utility.hpp"

using std::array;
using std::string;
using std::vector;
using base::utility::format;


SampleMass_FileBackingStore::SampleMass_FileBackingStore(string baseName)
    : FileBackingStore(baseName + ".massSamples"), membership(baseName + ".membership")
{
    if(!membership)
    {
        throw std::runtime_error(baseName + ".membership was not available for writing.");
    }
}

void SampleMass_FileBackingStore::save(vector<SampleMassRecord> data)
{
    if (!data.empty() && data.at(0).iter.val == 1)
    {
        header(data);
    }

    for ( auto d : data )
    {
        fout << base::utility::format << d.primaryMass
             << base::utility::format << d.massRatio;

        membership << base::utility::format << d.clusterMembership;
    }

    fout       << endl;
    membership << endl;
}


void SampleMass_FileBackingStore::header(vector<SampleMassRecord> data)
{
    for ( auto d : data )
    {
        fout << base::utility::format << d.starId + " mass"
             << base::utility::format << d.starId + " ratio";
    }

    fout << endl;
}
