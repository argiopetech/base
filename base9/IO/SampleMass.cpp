#include <array>
#include <string>

#include "SampleMass.hpp"
#include "Utility.hpp"

using std::array;
using std::string;
using std::vector;
using base::utility::format;


SampleMass_FileBackingStore::SampleMass_FileBackingStore(string baseName)
    : FileBackingStore(baseName + ".sampleMass.out")
{ ; }

void SampleMass_FileBackingStore::save(vector<SampleMassRecord> data)
{
    if (!data.empty() && data.front().iter.val == 1)
    {
        header(data);
    }

    for ( auto d : data )
    {
        fout << base::utility::format << d.iter.val
             << base::utility::format << d.starId
             << base::utility::format << d.primaryMass
             << base::utility::format << d.massRatio
             << base::utility::format << d.clusterMembership
             << endl;
    }
}


void SampleMass_FileBackingStore::header(vector<SampleMassRecord>)
{
    fout << base::utility::format << " iteration"
         << base::utility::format << " starId"
         << base::utility::format << " mass"
         << base::utility::format << " massRatio"
         << base::utility::format << " membership"
         << endl;
}
