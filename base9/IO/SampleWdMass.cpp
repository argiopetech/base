#include <array>
#include <string>

#include "SampleWdMass.hpp"
#include "Utility.hpp"

using std::array;
using std::string;
using std::vector;
using base::utility::format;


SampleWdMass_FileBackingStore::SampleWdMass_FileBackingStore(string baseName)
    : FileBackingStore(baseName + ".wd.mass"), membership(baseName + ".wd.membership"),
      precLogAge(baseName + ".wd.precLogAge"), coolingAge(baseName + ".wd.coolingAge"),
      logTeff(baseName + ".wd.logTeff"), logG(baseName + ".wd.logg")
{
    if(!membership)
    {
        throw std::runtime_error(baseName + ".wd.membership was not available for writing.");
    }

    if(!precLogAge)
    {
        throw std::runtime_error(baseName + ".wd.precLogAge was not available for writing.");
    }

    if(!coolingAge)
    {
        throw std::runtime_error(baseName + ".wd.coolingAge was not available for writing.");
    }

    if(!logTeff)
    {
        throw std::runtime_error(baseName + ".wd.logTeff was not available for writing.");
    }

    if(!logG)
    {
        throw std::runtime_error(baseName + ".wd.logg was not available for writing.");
    }
}

void SampleWdMass_FileBackingStore::save(vector<SampleWdMassRecord> data)
{
    for ( auto d : data )
    {
            fout       << base::utility::format << d.mass;
            membership << base::utility::format << d.clusterMembership;
            precLogAge << base::utility::format << d.precursorLogAge;
            coolingAge << base::utility::format << d.coolingAge;
            logTeff    << base::utility::format << d.logTeff;
            logG       << base::utility::format << d.logG;
    }

    fout       << endl;
    membership << endl;
    precLogAge << endl;
    coolingAge << endl;
    logTeff    << endl;
    logG       << endl;
}


void SampleWdMass_FileBackingStore::header(SampleWdMassRecord)
{ ; }
