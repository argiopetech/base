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
