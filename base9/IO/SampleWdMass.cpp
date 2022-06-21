#include <array>
#include <string>

#include "SampleWdMass.hpp"
#include "Utility.hpp"

using std::array;
using std::string;
using std::vector;
using base::utility::format;


SampleWdMass_FileBackingStore::SampleWdMass_FileBackingStore(string baseName)
    : FileBackingStore(baseName + ".sampleWDMass.out")
{ ; }

void SampleWdMass_FileBackingStore::save(vector<SampleWdMassRecord> data)
{
    if (!data.empty() && data.front().iter.val == 1)
    {
        header(data);
    }

    for ( auto d : data )
    {
        fout << base::utility::format          << d.iter.val          << ' '
             << std::setw(longestStarIdLength) << d.starId            << ' '
             << base::utility::format          << d.mass              << ' '
             << base::utility::format          << d.clusterMembership << ' '
             << base::utility::format          << d.precursorLogAge   << ' '
             << base::utility::format          << d.coolingAge        << ' '
             << base::utility::format          << d.logTeff           << ' '
             << base::utility::format          << d.logLittleG
             << endl;
    }
}


void SampleWdMass_FileBackingStore::header(std::vector<SampleWdMassRecord> records)
{
    for (auto record : records)
    {
        longestStarIdLength = std::max(longestStarIdLength, record.starId.size());
    }

     fout << base::utility::format          << "iteration"  << ' '
          << std::setw(longestStarIdLength) << "starId"     << ' '
          << base::utility::format          << "mass"       << ' '
          << base::utility::format          << "clustMembr" << ' '
          << base::utility::format          << "precLogAge" << ' '
          << base::utility::format          << "coolLogAge" << ' '
          << base::utility::format          << "logTeff"    << ' '
          << base::utility::format          << "logg"
          << endl;
}
