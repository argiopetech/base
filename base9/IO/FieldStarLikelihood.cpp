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
