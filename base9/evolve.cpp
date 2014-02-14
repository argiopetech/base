#include <array>
#include <vector>

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "Cluster.hpp"
#include "Star.hpp"

#include "evolve.hpp"
#include "marg.hpp"
#include "MsFilterSet.hpp"
#include "WhiteDwarf.hpp"

using std::array;
using std::vector;

/**************************************************************************************
 * This routine is the master subroutine of the simulation code and is organized so mcmc.c
 * can call this routine one pair of stars at a time for any of a range of hypothetical
 * stellar, cluster, or model properties.  This routine in turn calls other subroutines.
 * The parameters used by this routine are the cluster properties of age, metallicity, distance,
 * reddening, the model set (which indicates a combination of a stellar evolution model,
 * an initial-final mass relation, a WD cooling model, and a WD atmosphere model),  and the
 * ZAMS masses of the two stars.  This routine does not control binary fraction (parent
 * routines will take care of binaries by creating a primary and a secondary mass, the
 * latter of which can be zero), nor whether a WD is a DA or a DB (again, controlled by
 * parent routine).  Using all of these inputs, this routine updates the photometry of the
 * star structure in the whichever of the U through K filters are being used (stored in useFilt).
***************************************************************************************/
void evolve (const Cluster &pCluster, const Model &evoModels, array<double, FILTS> &globalMags, const vector<int> &filters,  Star &star, array<double, 2> &ltau)
{
    array<double, 2> mass;
    Matrix<double, 3, FILTS> mag;

    // AGBt_zmass never set because age and/or metallicity out of range of models.
    if (pCluster.AGBt_zmass < EPS)
    {
        throw WDBoundsError("Bounds error in evolve.cpp");
    }

    mass[0] = star.getMass1();
    mass[1] = star.getMass2();

    for (int cmpnt = 0; cmpnt < 2; cmpnt++)
    {
        mag[cmpnt] = star.setMags(cmpnt, mass[cmpnt], pCluster, evoModels, filters, ltau[cmpnt]);
    }

    star.deriveCombinedMags(mag, pCluster, evoModels, filters);
}
