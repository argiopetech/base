#include <array>
#include <string>

#include "SampleWdMass.hpp"
#include "Utility.hpp"

using std::array;
using std::string;
using std::vector;
using base::utility::format;


SampleWdMass_FileBackingStore::SampleWdMass_FileBackingStore(string baseName)
    : FileBackingStore(baseName + ".sampleWdMass.out")
{ ; }

void SampleWdMass_FileBackingStore::save(vector<SampleWdMassRecord> data)
{
    if (!data.empty() && data.front().iter.val == 1)
    {
        header(data.front());
    }

    for ( auto d : data )
    {
        fout << base::utility::format << d.iter.val
             << base::utility::format << d.starId
             << base::utility::format << d.mass
             << base::utility::format << d.clusterMembership
             << base::utility::format << d.precursorLogAge
             << base::utility::format << d.coolingAge
             << base::utility::format << d.logTeff
             << base::utility::format << d.logG
             << endl;
    }
}


void SampleWdMass_FileBackingStore::header(SampleWdMassRecord)
{
    const array<string, 8> paramNames = { "  iteration",
                                          "     starId",
                                          "       mass",
                                          " clustMembr",
                                          " precLogAge",
                                          " coolLogAge",
                                          "    logTeff",
                                          "       logG"};

    for (auto p : paramNames)
    {
        fout << p;
    }

    fout << endl;
}
