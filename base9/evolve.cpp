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
array<double, FILTS> evolve (const Cluster &pCluster, const Model &evoModels, const vector<int> &filters,  StellarSystem &system)
{
    Matrix<double, 2, FILTS> mag;

    // AGBt_zmass never set because age and/or metallicity out of range of models.
    if (pCluster.AGBt_zmass < EPS)
    {
        throw WDBoundsError("Bounds error in evolve.cpp");
    }

    mag.at(0) = system.primary.getMags(system.primary.mass, pCluster, evoModels, filters);
    mag.at(1) = system.secondary.getMags(system.secondary.mass, pCluster, evoModels, filters);

    return deriveCombinedMags(mag.at(0), mag.at(1), pCluster, evoModels, filters);
}
