#include <array>

#include "constants.hpp"
#include "FilterSet.hpp"

using std::array;

UBVRIJHK::UBVRIJHK()
    : FilterSet({"U", "B", "V", "R", "I", "J", "H", "K"})
{}

array<double, FILTS> UBVRIJHK::calcAbsCoeffs() const
{
    return { 1.569  // Cardelli, Clayton, Mathis 1989, table 3
           , 1.337  // yields A_u -> A_k = f(A_v), for standard filters
           , 1.0
           , 0.751
           , 0.479
           , 0.282
           , 0.190
           , 0.114};
}


SDSS::SDSS()
    : FilterSet({"u", "g", "r", "i", "z", "J", "H", "K"})
{}

array<double, FILTS> SDSS::calcAbsCoeffs() const
{
    return { 5.155 / 3.1    // Stoughton et al. (2002, AJ, 123, 485)
           , 3.793 / 3.1    // Table 22, which gives Afilter/E(B-V )
           , 2.751 / 3.1    // We use Afilter/Av, so all are divided
           , 2.086 / 3.1    // by Rv = Av/E(B-V) = 3.1, consistent
           , 1.479 / 3.1    // Cardelli et al. (1989).
           , 0.282  // JHK come from Cardelli (see above)
           , 0.190
           , 0.114};
}


ACS::ACS()
    : FilterSet({"F435W","F475W","F550M","F555W","F606W","F625W","F775W","F814W"})
{}

array<double, FILTS> ACS::calcAbsCoeffs() const
{
    return { 4.081 / 3.1    // they used R=3.1; also derived via Cardelli et al. values
           , 3.634 / 3.1    // Ext. Ratios A(P)/E(B-V) in ACS/WFC Filters for diff. SEDs
           , 3.042 / 3.1    // SED F435W F475W F550M F555W F606W F625W F775W F814W
           , 3.177 / 3.1    // O5  4.192 3.773 3.052 3.233 2.936 2.673 2.005 1.864
           , 2.809 / 3.1    // G2  4.081 3.634 3.042 3.177 2.809 2.637 1.982 1.825
           , 2.637 / 3.1    // M0  3.994 3.555 3.030 3.115 2.716 2.616 1.965 1.796
           , 1.982 / 3.1
           , 1.825 / 3.1};    // using values for G2 star, usually ~2% of O5 or M0 value
}


UVIS::UVIS()
    : FilterSet({"UVf275w", "UVf336w", "UVf438w"})
{}

array<double, FILTS> UVIS::calcAbsCoeffs() const
{
    return { 2.094     // Based on calculations of Ax/Av by Rachel Wagner-Kaiser
           , 1.665     // from the theoretical values given in the Cardelli paper.
           , 1.343};
}
