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
        fout << base::utility::format          << d.iter.val          << ' '
             << std::setw(longestStarIdLength) << d.starId            << ' '
             << base::utility::format          << d.primaryMass       << ' '
             << base::utility::format          << d.massRatio         << ' '
             << base::utility::format          << d.clusterMembership
             << endl;
    }
}


void SampleMass_FileBackingStore::header(vector<SampleMassRecord> records)
{
    for (auto record : records)
    {
        longestStarIdLength = std::max(longestStarIdLength, record.starId.size());
    }

    fout << base::utility::format          << "iteration"  << ' '
         << std::setw(longestStarIdLength) << "starId"     << ' '
         << base::utility::format          << "mass"       << ' '
         << base::utility::format          << "massRatio"  << ' '
         << base::utility::format          << "membership"
         << endl;
}
