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
