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

void SinglePopMcmc_FileBackingStore::save(SinglePopMcmcRecord data)
{
    const auto &clust = data.clust;

    if (data.iter.val == 1)
    {
        header(clust.priorVar, data.modIsParallax);
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


void SinglePopMcmc_FileBackingStore::header(array<double, NPARAMS> const &priors, bool modIsParallax)
{
    auto distanceName = modIsParallax ? "    parallax" : "     modulus";

    const array<string, NPARAMS> paramNames = {  "     logAge",
                                                "           Y",
                                                "         FeH",
                                                  distanceName,
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
