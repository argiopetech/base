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


SampleMass_SqlBackingStore::SampleMass_SqlBackingStore(const RunData &bootstrap)
    : SqlBackingStore(bootstrap)
{ buildInsertStatement(); }


SampleMass_SqlBackingStore::SampleMass_SqlBackingStore(std::string filename)
    : SqlBackingStore(filename)
{ buildInsertStatement(); }

SampleMass_SqlBackingStore::~SampleMass_SqlBackingStore()
{
    sqlite3_finalize(insert);
    insert = nullptr;
}

void SampleMass_SqlBackingStore::buildInsertStatement()
{
    dbPrepare( "insert into sample_mass values (?1, ?2, ?3, ?4, ?5, ?6, ?7);",
               &insert, "Preparing sample mass insert");
}

void SampleMass_SqlBackingStore::save(vector<SampleMassRecord> data)
{
    beginTransaction("sample mass");

    for (auto d : data)
    {
        sqlite3_bind_int(insert, 1, run);
        sqlite3_bind_int(insert, 2, d.referencedRun);
        sqlite3_bind_int(insert, 3, d.iter.val);
        sqlite3_bind_text(insert, 4, d.starId.data(), d.starId.length(), SQLITE_STATIC);

        sqlite3_bind_double(insert, 5, d.primaryMass);
        sqlite3_bind_double(insert, 6, d.massRatio);
        sqlite3_bind_double(insert, 7, d.clusterMembership);

        dbStepAndReset(insert, "Inserting record into sample_mass");
    }

    endTransaction("sample mass");
}
