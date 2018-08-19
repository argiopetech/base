#include <array>
#include <string>

#include "FieldStarLikelihood.hpp"
#include "Utility.hpp"

using std::array;
using std::string;
using base::utility::format;


FieldStarLikelihood_FileBackingStore::FieldStarLikelihood_FileBackingStore(string baseName)
    : FileBackingStore(baseName + ".fslikelihood")
{ ; }

void FieldStarLikelihood_FileBackingStore::save(FieldStarLikelihoodRecord data)
{
    header(data.fsLike);
}


void FieldStarLikelihood_FileBackingStore::header(double fsLike)
{
    fout << format << fsLike << std::endl;
}


FieldStarLikelihood_SqlBackingStore::FieldStarLikelihood_SqlBackingStore(const RunData &bootstrap)
    : SqlBackingStore(bootstrap)
{ buildInsertStatement(); }



FieldStarLikelihood_SqlBackingStore::FieldStarLikelihood_SqlBackingStore(std::string filename)
    : SqlBackingStore(filename)
{ buildInsertStatement(); }


FieldStarLikelihood_SqlBackingStore::~FieldStarLikelihood_SqlBackingStore()
{
    sqlite3_finalize(insert);
    insert = nullptr;

}


void FieldStarLikelihood_SqlBackingStore::buildInsertStatement()
{
    dbPrepare( "insert into field_star_likelihood values (?1, ?2);",
               &insert, "Preparing field_star_likelihood insert");
}


void FieldStarLikelihood_SqlBackingStore::save(FieldStarLikelihoodRecord data)
{
    sqlite3_bind_int   (insert, 1, run);
    sqlite3_bind_double(insert, 2, data.fsLike);

    dbStepAndReset(insert, "Inserting field star likelihood");
}
