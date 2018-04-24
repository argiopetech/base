#include <array>
#include <string>

#include "SinglePopMcmc.hpp"
#include "Utility.hpp"

using std::array;
using std::string;
using base::utility::format;


SinglePopMcmc_FileBackingStore::SinglePopMcmc_FileBackingStore(string baseName)
    : FileBackingStore(baseName + ".res")
{ ; }

void SinglePopMcmc_FileBackingStore::save(Iteration iter, SinglePopMcmcRecord data)
{
    const auto &clust = data.clust;

    if (iter.val == 1)
    {
        header(clust.priorVar);
    }

    for (int p = 0; p < NPARAMS; p++)
    {
        if (clust.priorVar.at(p) > EPS || p == FEH || p == MOD || p == ABS)
        {
            fout << format << clust.getParam(p) << ' ';
        }
    }

    fout << format << data.logPost
         << format << static_cast<int>(data.stage) << endl;
}


void SinglePopMcmc_FileBackingStore::header(array<double, NPARAMS> const &priors)
{
    const array<string, NPARAMS> paramNames = {  "     logAge",
                                                "           Y",
                                                "         FeH",
                                                "     modulus",
                                                "  absorption",
                                                " carbonicity",
                                                "   IFMRconst",
                                                "     IFMRlin",
                                                "    IFMRquad"};

    for (int p = 0; p < NPARAMS; p++)
    {
        if (priors.at(p) > EPSILON || p == MOD || p == FEH || p == ABS)
        {
            fout << paramNames.at(p);
        }
    }

    fout << "      logPost" << "     stage" << endl;
}


SinglePopMcmc_SqlBackingStore::SinglePopMcmc_SqlBackingStore(const RunData &bootstrap)
    : SqlBackingStore(bootstrap)
{ buildInsertStatement(); }


SinglePopMcmc_SqlBackingStore::SinglePopMcmc_SqlBackingStore(std::string filename)
    : SqlBackingStore(filename)
{ buildInsertStatement(); }

SinglePopMcmc_SqlBackingStore::~SinglePopMcmc_SqlBackingStore()
{
    sqlite3_finalize(insert);
    insert = nullptr;
}

void SinglePopMcmc_SqlBackingStore::buildInsertStatement()
{
    dbPrepare( "insert into single_pop values (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13);",
               &insert, "Preparing single pop insert");
}

void SinglePopMcmc_SqlBackingStore::save(Iteration iter, SinglePopMcmcRecord data)
{
    auto &clust = data.clust;

    sqlite3_bind_int(insert, 1, run);
    sqlite3_bind_int(insert, 2, iter.val);

    (clust.priorVar[AGE] > EPS)
        ? sqlite3_bind_double(insert, 3, clust.age)
        : sqlite3_bind_null(insert, 3);

    (clust.priorVar[YYY] > EPS)
        ? sqlite3_bind_double(insert, 4, clust.yyy)
        : sqlite3_bind_null(insert, 4);

    sqlite3_bind_double(insert, 5, clust.feh);
    sqlite3_bind_double(insert, 6, clust.mod);
    sqlite3_bind_double(insert, 7, clust.abs);

    (clust.priorVar[CARBONICITY] > EPS)
        ? sqlite3_bind_double(insert, 8, clust.carbonicity)
        : sqlite3_bind_null(insert, 8);

    (clust.priorVar[IFMR_INTERCEPT] > EPS)
        ? sqlite3_bind_double(insert, 9, clust.ifmrIntercept)
        : sqlite3_bind_null(insert, 9);

    (clust.priorVar[IFMR_SLOPE] > EPS)
        ? sqlite3_bind_double(insert, 10, clust.ifmrSlope)
        : sqlite3_bind_null(insert, 10);

    (clust.priorVar[IFMR_QUADCOEF] > EPS)
        ? sqlite3_bind_double(insert, 11, clust.ifmrQuadCoef)
        : sqlite3_bind_null(insert, 11);

    sqlite3_bind_double(insert, 12, data.logPost);
    sqlite3_bind_int   (insert, 13, static_cast<int>(data.stage));

    dbStepAndReset(insert, "Inserting record into single_pop");
}
