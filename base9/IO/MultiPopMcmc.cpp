#include <array>
#include <string>

#include "MultiPopMcmc.hpp"
#include "Utility.hpp"

using std::array;
using std::string;
using base::utility::format;

MultiPopMcmc_FileBackingStore::MultiPopMcmc_FileBackingStore(string baseName)
    : FileBackingStore(baseName + ".res")
{ ; }

void MultiPopMcmc_FileBackingStore::save(MultiPopMcmcRecord data)
{
    const auto &clust0 = data.clust0;
    const auto &clust1 = data.clust1;

    if (data.iter.val == 1)
    {
        header(clust0.priorVar, data.modIsParallax);
    }

    for (int p = 0; p < NPARAMS; p++)
    {
        if (clust0.priorVar.at(p) > EPS || p == FEH || p == MOD || p == ABS)
        {
            fout << format << clust0.getParam(p) << ' ';
        }
    }

    fout << format << clust0.yyy << ' '
         << format << clust1.yyy << ' '
         << format << data.lambda << ' ';

    fout << format << data.logPost
         << format << static_cast<int>(data.stage) << endl;
}


void MultiPopMcmc_FileBackingStore::header(array<double, NPARAMS> const &priors, bool modIsParallax)
{
    auto distanceName = modIsParallax ? "    parallax" : "     modulus";

    const array<string, NPARAMS> paramNames = {  "     logAge",
                                                "          YA",
                                                "         FeH",
                                                  distanceName,
                                                "  absorption",
                                                " carbonicity",
                                                "          YB",
                                                "      Lambda",
                                                "    IFMRquad"};

    for (int p = 0; p < NPARAMS; p++)
    {
        if ((priors.at(p) > EPSILON || p == MOD || p == FEH || p == ABS) && !(p == YYA || p == YYB || p == LAMBDA))
        {
            fout << paramNames.at(p);
        }
    }
    fout << paramNames.at(YYA) << paramNames.at(YYB) << paramNames.at(LAMBDA);
    fout << "      logPost" << "     stage" << endl;
}


MultiPopMcmc_SqlBackingStore::MultiPopMcmc_SqlBackingStore(const RunData &bootstrap)
    : SqlBackingStore(bootstrap)
{ buildInsertStatement(); }


MultiPopMcmc_SqlBackingStore::MultiPopMcmc_SqlBackingStore(std::string filename)
    : SqlBackingStore(filename)
{ buildInsertStatement(); }

MultiPopMcmc_SqlBackingStore::~MultiPopMcmc_SqlBackingStore()
{
    sqlite3_finalize(insert);
    insert = nullptr;
}

void MultiPopMcmc_SqlBackingStore::buildInsertStatement()
{
    dbPrepare( "insert into multi_pop values (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12);",
               &insert, "Preparing multi pop insert");
}

void MultiPopMcmc_SqlBackingStore::save(MultiPopMcmcRecord data)
{
    const auto &clust0 = data.clust0;
    const auto &clust1 = data.clust1;

    sqlite3_bind_int(insert, 1, run);
    sqlite3_bind_int(insert, 2, data.iter.val);

    (clust0.priorVar[AGE] > EPS)
        ? sqlite3_bind_double(insert, 3, clust0.age)
        : sqlite3_bind_null(insert, 3);

    sqlite3_bind_double(insert, 4, clust0.feh);
    sqlite3_bind_double(insert, 5, clust0.mod);
    sqlite3_bind_double(insert, 6, clust0.abs);

    (clust0.priorVar[CARBONICITY] > EPS)
        ? sqlite3_bind_double(insert, 7, clust0.carbonicity)
        : sqlite3_bind_null(insert, 7);

    sqlite3_bind_double(insert,  8, clust0.yyy);
    sqlite3_bind_double(insert,  9, clust1.yyy); // Currently the only piece of duplicated data
    sqlite3_bind_double(insert, 10, data.lambda);
    sqlite3_bind_double(insert, 11, data.logPost);
    sqlite3_bind_int   (insert, 12, static_cast<int>(data.stage));

    dbStepAndReset(insert, "Inserting record into multi_pop");
}
