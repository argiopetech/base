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

void MultiPopMcmc_FileBackingStore::save(Iteration iter, MultiPopMcmcRecord data)
{
    const auto &clust0 = data.clust0;
    const auto &clust1 = data.clust1;
    
    if (iter.val == 1)
    {
        header(clust0.priorVar);
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


void MultiPopMcmc_FileBackingStore::header(array<double, NPARAMS> const &priors)
{
    const array<string, NPARAMS> paramNames = {  "     logAge",
                                                "          YA",
                                                "         FeH",
                                                "     modulus",
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
