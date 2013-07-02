/*** Decides whether to accept proposed Metropolis-Hastings steps ***/

#include <vector>

#include <cmath>
#include <cstdio>

#include "samplers.hpp"
#include "densities.hpp"
#include "decide.hpp"
#include "mt19937ar.hpp"
#include "evolve.hpp"
#include "structures.hpp"

using std::vector;

/*** Decides whether to accept a proposed jump between field star and cluster star models ***/
void decideFieldStar (Star stars1[], Cluster &pCluster, FILE * wFile, Model &evoModels)
{
    int j;
    double u, alpha, post1, post2;
    vector<Star> stars2(pCluster.nStars);

    for (j = 0; j < pCluster.nStars; j++)
    {                           // For each star,
        stars2.at(j) = stars1[j];  // copy stars to new array,
        propFieldStar (&(stars2.at(j)));   // propose a new field star status,

        post1 = logPost1Star (stars1[j], pCluster, evoModels);
        post2 = logPost1Star (stars2.at(j), pCluster, evoModels);

        if (fabs (post2 + HUGE_VAL) < EPS)
            continue;                   // Proposed star no good.  Leave stars1 alone.

        // Don't bother with this if the prior and proposal densities will cancel out
        if (fabs (stars1[j].clustStarPriorDens - stars1[j].clustStarProposalDens) > EPS)
        {
            // Multiply by prior probability of star's FS status
            if (stars1[j].isFieldStar)
                post1 += log (1 - stars1[j].clustStarPriorDens);
            else
                post1 += log (stars1[j].clustStarPriorDens);
            if (stars2.at(j).isFieldStar)
                post2 += log (1 - stars2.at(j).clustStarPriorDens);
            else
                post2 += log (stars2.at(j).clustStarPriorDens);
            // Divide by proposal densities
            if (stars1[j].isFieldStar)
                post1 -= log (1 - stars1[j].clustStarProposalDens);
            else
                post1 -= log (stars1[j].clustStarProposalDens);
            if (stars2.at(j).isFieldStar)
                post2 -= log (1 - stars2.at(j).clustStarProposalDens);
            else
                post2 -= log (stars2.at(j).clustStarProposalDens);
        }

        alpha = post2 - post1;
        u = genrand_res53 ();
        if (u < 1.e-15)
            u = 1.e-15;
        u = log (u);

        if (u < alpha)
            stars1[j].isFieldStar = stars2.at(j).isFieldStar;      // Accept proposed star (copy back to the stars1 array)
    }
}

/*** Decides whether to accept a proposed mass ***/
void decideMass (Chain &mc, Model &evoModels)
{
    int j;
    double u, alpha, post1, post2;
    vector<Star> stars2(mc.clust.nStars);

    for (j = 0; j < mc.clust.nStars; j++)
    {                           // For each star,
        stars2.at(j) = mc.stars.at(j);       // copy stars to new array,
        propMass (&stars2.at(j));  // propose a new mass,
        stars2.at(j).boundsFlag = 0;       // and set the boundsFlag to zero
    }

    for (auto s : stars2)
        evolve (mc.clust, evoModels, s);    // Evolve all the (proposed) stars at once

    for (j = 0; j < mc.clust.nStars; j++)
    {                           // Accept or reject each star individually
        if (stars2.at(j).boundsFlag || getMass1 (stars2.at(j), mc.clust) < EPS)
            mc.rejectMass[j]++;        // Proposed star no good.  Leave stars1 alone.
        else
        {
            post2 = logPost1Star (stars2.at(j), mc.clust, evoModels);
            if (fabs (post2 + HUGE_VAL) < EPS)
                mc.rejectMass[j]++;    // Proposed star no good.  Leave stars1 alone.
            else
            {
                post1 = logPost1Star (mc.stars.at(j), mc.clust, evoModels);
                alpha = post2;
                alpha -= post1;

                u = genrand_res53 ();
                if (u < 1.e-15)
                    u = 1.e-15;
                u = log (u);

                if (u < alpha)
                {
                    mc.acceptMass[j]++;        // Accept proposed star
                    mc.stars.at(j) = stars2.at(j);   // And copy back to the stars1 array
                }
                else
                    mc.rejectMass[j]++;        // Proposed star no good.  Leave stars1 alone.
            }
        }
    }
}

/*** Decides whether to accept a proposed mass ratio ***/
void decideMassRatio (Chain &mc, Model &evoModels)
{
    double u, alpha, post1, post2;
    int j;

    vector<Star> stars2(mc.clust.nStars);

    for (j = 0; j < mc.clust.nStars; j++)
    {                           // For each star,
        stars2.at(j) = mc.stars.at(j);       // copy stars to new array,
        propMassRatio (&stars2.at(j));     // propose a new mass ratio,
        stars2.at(j).boundsFlag = 0;       // and set the boundsFlag to zero
    }

    for (auto s : stars2)
        evolve (mc.clust, evoModels, s);    // Evolve all the (proposed) stars at once

    for (j = 0; j < mc.clust.nStars; j++)
    {                           // Accept or reject each star individually
        if (stars2.at(j).boundsFlag || getMass1 (stars2.at(j), mc.clust) < EPS)
            mc.rejectMassRatio[j]++;   // Proposed star no good.  Leave stars1 alone.
        else if (stars2.at(j).massRatio > 1.0 || stars2.at(j).massRatio < 0.0)
            mc.rejectMassRatio[j]++;   // Proposed star no good.  Leave stars1 alone.
        else
        {
            post2 = logPost1Star (stars2.at(j), mc.clust, evoModels);
            if (fabs (post2 + HUGE_VAL) < EPS)
                mc.rejectMassRatio[j]++;       // Proposed star no good.  Leave stars1 alone.
            else
            {
                post1 = logPost1Star (mc.stars.at(j), mc.clust, evoModels);
                alpha = post2;
                alpha -= post1;

                u = genrand_res53 ();
                if (u < 1.e-15)
                    u = 1.e-15;
                u = log (u);

                if (u < alpha)
                {
                    mc.acceptMassRatio[j]++;   // Accept proposed star
                    mc.stars.at(j) = stars2.at(j);   // And copy back to the stars1 array
                }
                else
                    mc.rejectMassRatio[j]++;   // Proposed star no good.  Leave stars1 alone.
            }
        }
    }
}

// Decides whether to accept a proposed cluster property
Cluster decideClust (Cluster clust1, Star stars1[], const int FS_ON_STATE, int *accept, int *reject, const int SAMPLE_TYPE, FILE * w_ptr, Model &evoModels)
{
    int j;
    double u, alpha, post1 = 0.0, post2 = 0.0;

    vector<Star> stars2(clust1.nStars);
    Cluster clust2;

    clust2 = clust1;
    propClustParam (&clust2, SAMPLE_TYPE);      // propose a new value

    post1 = logPriorClust (clust1, evoModels);
    post2 = logPriorClust (clust2, evoModels);

    if (fabs (post2 + HUGE_VAL) < EPS)
    {
        (*reject)++;
        return clust1;
    }

    for (j = 0; j < clust1.nStars; j++)
    {
        if (getMass1 (stars1[j], clust2) < EPS)
        {
            (*reject)++;
            return clust1;
        }
        stars2.at(j) = stars1[j];
        stars2.at(j).boundsFlag = 0;
    }

    for( auto s : stars2)
        evolve (clust2, evoModels, s);

    for (j = 0; j < clust1.nStars; j++)
    {
        if (stars2.at(j).boundsFlag)
        {
            (*reject)++;
            return clust1;
        }
        if (FS_ON_STATE || stars1[j].useDuringBurnIn)
        {

            post1 += logPost1Star (stars1[j], clust1, evoModels);
            post2 += logPost1Star (stars2.at(j), clust2, evoModels);

            if (fabs (post2 + HUGE_VAL) < EPS)
            {
                (*reject)++;
                return clust1;
            }
        }
    }
    alpha = (post2 - post1);

    u = genrand_res53 ();
    if (u < 1.e-15)
        u = 1.e-15;
    u = log (u);

    if (u < alpha)
    {
        for (j = 0; j < clust1.nStars; j++)
            stars1[j] = stars2.at(j);
        (*accept)++;
        return clust2;
    }
    else
    {
        (*reject)++;
        return clust1;
    }
}

// Draw a new varScale from a scaled Inv-gamma distribution.
// This assumes that the prior distribution for the varScale parameter is Inv-chisq(prior_df)
void updateVarScale (Star stars[], Cluster &pCluster, Model &evoModels)
{
    int nClustStars = 0, i, j;
    double scale = 0.0;
    double prior_df = 3.0;

    for (i = 0; i < pCluster.nStars; i++)
    {
        if (!stars[i].isFieldStar)
        {
            nClustStars++;
            for (j = 0; j < evoModels.numFilts; j++)
            {
                if (stars[i].variance[j] > 0)
                    scale += 0.5 * ((stars[i].photometry[j] - stars[i].obsPhot[j]) * (stars[i].photometry[j] - stars[i].obsPhot[j]) / stars[i].variance[j]);
            }
        }
    }
    scale += 0.5;
    double g = gamdev ((evoModels.numFilts * nClustStars + prior_df) / 2.0);

    pCluster.varScale = scale / g;
}
