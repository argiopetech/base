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


SampleWdMass_SqlBackingStore::SampleWdMass_SqlBackingStore(const RunData &bootstrap)
    : SqlBackingStore(bootstrap)
{ buildInsertStatement(); }


SampleWdMass_SqlBackingStore::SampleWdMass_SqlBackingStore(std::string filename)
    : SqlBackingStore(filename)
{ buildInsertStatement(); }

SampleWdMass_SqlBackingStore::~SampleWdMass_SqlBackingStore()
{
    sqlite3_finalize(insert);
    insert = nullptr;
}

void SampleWdMass_SqlBackingStore::buildInsertStatement()
{
    dbPrepare( "insert into sample_wd_mass values (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10);",
               &insert, "Preparing sample WD mass insert");
}

void SampleWdMass_SqlBackingStore::save(vector<SampleWdMassRecord> data)
{
    beginTransaction("sample WD mass");

    for (auto d : data)
    {
        sqlite3_bind_int(insert, 1, run);
        sqlite3_bind_int(insert, 2, d.referencedRun);
        sqlite3_bind_int(insert, 3, d.iter.val);
        sqlite3_bind_text(insert, 4, d.starId.data(), d.starId.length(), SQLITE_STATIC);

        sqlite3_bind_double(insert, 5, d.mass);
        sqlite3_bind_double(insert, 6, d.clusterMembership);
        sqlite3_bind_double(insert, 7, d.precursorLogAge);
        sqlite3_bind_double(insert, 8, d.coolingAge);
        sqlite3_bind_double(insert, 9, d.logTeff);
        sqlite3_bind_double(insert, 10, d.logG);

        dbStepAndReset(insert, "Inserting record into sample_wd_mass");
    }

    endTransaction("sample WD mass");
}
