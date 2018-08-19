#include <array>
#include <string>

#include "Star.hpp"

using std::string;


Star_SqlBackingStore::Star_SqlBackingStore(const RunData &bootstrap)
    : SqlBackingStore(bootstrap)
{ buildInsertStatement(); }


Star_SqlBackingStore::Star_SqlBackingStore(std::string filename)
    : SqlBackingStore(filename)
{ buildInsertStatement(); }

Star_SqlBackingStore::~Star_SqlBackingStore()
{
    sqlite3_finalize(insertStar);
    insertStar = nullptr;

    sqlite3_finalize(insertPhotometry);
    insertPhotometry = nullptr;
}

void Star_SqlBackingStore::buildInsertStatement()
{
    dbPrepare( "insert into star values (?1, ?2, ?3, ?4, ?5, ?6, ?7);",
               &insertStar, "Preparing star insert");

    dbPrepare( "insert into photometry values (?1, ?2, ?3, ?4, ?5);",
               &insertPhotometry, "Preparing photometry insert");
}

void Star_SqlBackingStore::save(StarRecord data)
{
    beginTransaction("photometry");

    {
        sqlite3_bind_int(insertStar, 1, run);
        sqlite3_bind_text(insertStar, 2, data.starId.data(), data.starId.length(), SQLITE_STATIC);
        sqlite3_bind_double(insertStar, 3, data.primaryMass);
        sqlite3_bind_double(insertStar, 4, data.secondaryMassRatio);
        sqlite3_bind_int(insertStar, 5, data.stage);
        sqlite3_bind_double(insertStar, 6, data.clusterMembershipPrior);
        sqlite3_bind_int(insertStar, 7, data.useDuringBurnin);

        dbStepAndReset(insertStar, "Inserting record into star");
    }

    for (auto phot : data.photometryRecords)
    {
        sqlite3_bind_int(insertPhotometry, 1, run);
        sqlite3_bind_text(insertPhotometry, 2, data.starId.data(), data.starId.length(), SQLITE_STATIC);

        sqlite3_bind_text(insertPhotometry, 3, phot.filter.data(), phot.filter.length(), SQLITE_STATIC);
        sqlite3_bind_double(insertPhotometry, 4, phot.magnitude);
        sqlite3_bind_double(insertPhotometry, 5, phot.stdDev);

        dbStepAndReset(insertPhotometry, "Inserting record into photometry");
    }

    endTransaction("photometry");
}
