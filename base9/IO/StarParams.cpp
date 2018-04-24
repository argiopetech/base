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

void StarParams_FileBackingStore::save(Iteration iter, StarParamsRecord data)
{
    const auto &starData = data.starData;

    if (iter.val == 1)
    {
        header(data.fsLike);
    }

    fout << format << starData.at(0);

    for (size_t i = 1; i < starData.size(); ++i)
    {
        fout << ' ' << format << starData.at(i);
    }

    fout << std::endl;
}


void StarParams_FileBackingStore::header(double fsLike)
{
    fout << format << fsLike << std::endl;
}


StarParams_SqlBackingStore::StarParams_SqlBackingStore(const RunData &bootstrap)
    : SqlBackingStore(bootstrap)
{ buildInsertStatement(); }



StarParams_SqlBackingStore::StarParams_SqlBackingStore(std::string filename)
    : SqlBackingStore(filename)
{ buildInsertStatement(); }


StarParams_SqlBackingStore::~StarParams_SqlBackingStore()
{
    sqlite3_finalize(insertFsLike);
    insertFsLike = nullptr;

    sqlite3_finalize(insertStarData);
    insertStarData = nullptr;
}


void StarParams_SqlBackingStore::buildInsertStatement()
{
    dbPrepare( "insert into field_star_likelihood values (?1, ?2);",
               &insertFsLike, "Preparing field_star_likelihood insert");

    dbPrepare( "insert into star_posterior values (?1, ?2, ?3, ?4, ?5);",
               &insertStarData, "Preparing star_posterior insert");
}


void StarParams_SqlBackingStore::save(Iteration iter, StarParamsRecord data)
{
    const auto &starData = data.starData;

    if (iter.val == 1)
    {
        sqlite3_bind_int   (insertFsLike, 1, run);
        sqlite3_bind_double(insertFsLike, 2, data.fsLike);

        dbStepAndReset(insertFsLike, "Inserting field star likelihood");
    }

    dbErrorIf(execOnly("BEGIN TRANSACTION"), "StarData begin transaction");

    int starCounter = 1;

    for (size_t i = 1; i < starData.size(); i += 2)
    {
        sqlite3_bind_int(insertStarData, 1, run);
        sqlite3_bind_int(insertStarData, 2, iter.val);
        sqlite3_bind_int(insertStarData, 3, starCounter++);

        sqlite3_bind_double(insertStarData, 4, starData.at(i - 1));
        sqlite3_bind_double(insertStarData, 5, starData.at(i));

        dbStepAndReset(insertStarData, "Inserting record into star_params");
    }

    dbErrorIf(execOnly("END TRANSACTION"), "StarData end transaction");
}
