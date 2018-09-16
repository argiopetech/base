#include <array>
#include <string>

#include "StarParams.hpp"
#include "Utility.hpp"

using std::array;
using std::string;
using base::utility::format;


StarParams_FileBackingStore::StarParams_FileBackingStore(string baseName)
    : FileBackingStore(baseName + ".starParams")
{ ; }

void StarParams_FileBackingStore::save(StarParamsRecord data)
{
    const auto &starData = data.starData;

    fout << format << starData.at(0);

    for (size_t i = 1; i < starData.size(); ++i)
    {
        fout << ' ' << format << starData.at(i);
    }

    fout << std::endl;
}


StarParams_SqlBackingStore::StarParams_SqlBackingStore(const RunData &bootstrap)
    : SqlBackingStore(bootstrap)
{ buildInsertStatement(); }


StarParams_SqlBackingStore::StarParams_SqlBackingStore(const SqlBackingStore &bootstrap)
    : SqlBackingStore(bootstrap)
{ buildInsertStatement(); }


StarParams_SqlBackingStore::StarParams_SqlBackingStore(std::string filename)
    : SqlBackingStore(filename)
{ buildInsertStatement(); }


StarParams_SqlBackingStore::~StarParams_SqlBackingStore()
{
    sqlite3_finalize(insert);
    insert = nullptr;
}


void StarParams_SqlBackingStore::buildInsertStatement()
{
    dbPrepare( "insert into star_posterior values (?1, ?2, ?3, ?4, ?5);",
               &insert, "Preparing star_posterior insert");
}


void StarParams_SqlBackingStore::save(StarParamsRecord data)
{
    const auto &starData = data.starData;

    beginTransaction("individual star posteriors");

    for (size_t i = 1, s = 0; i < starData.size(); i += 2, ++s)
    {
        sqlite3_bind_int(insert, 1, run);
        sqlite3_bind_int(insert, 2, data.iter.val);
        sqlite3_bind_text(insert, 3, data.starNames[s].data(), data.starNames[s].length(), SQLITE_STATIC);

        sqlite3_bind_double(insert, 4, starData.at(i - 1));
        sqlite3_bind_double(insert, 5, starData.at(i));

        dbStepAndReset(insert, "Inserting individual star posterior");
    }

    endTransaction("individual star posteriors");
}
