#ifndef IO_RECORDS_HPP
#define IO_RECORDS_HPP

#include "BackingStore.hpp"


struct MultiPopMcmcRecord
{
    Iteration iter;

    const AdaptiveMcmcStage stage;

    const double lambda;

    const Cluster clust0;
    const Cluster clust1;

    const double logPost;

    const bool modIsParallax;
};

typedef BackingStore<MultiPopMcmcRecord> MultiPopBackingStore;


struct SinglePopMcmcRecord
{
    Iteration iter;

    const AdaptiveMcmcStage stage;
    const Cluster clust;
    const double logPost;

    const bool modIsParallax;
};

typedef BackingStore<SinglePopMcmcRecord> SinglePopBackingStore;


struct StarParamsRecord
{
    Iteration iter;

    const std::vector<string> starNames;
    const std::vector<double> starData;
};

typedef BackingStore<StarParamsRecord> StarParamsBackingStore;


struct PhotometryRecord
{
    const std::string starId;
    const std::string filter;
    const double magnitude;
    const double stdDev;
};

typedef BackingStore<PhotometryRecord> PhotometryBackingStore;


struct StarRecord
{
    const string starId;
    const double primaryMass;
    const double secondaryMassRatio;
    const int stage;
    const double clusterMembershipPrior;
    const bool useDuringBurnin;
    const std::vector<PhotometryRecord> photometryRecords;
};

typedef BackingStore<StarRecord> StarBackingStore;

struct FieldStarLikelihoodRecord
{
    const double fsLike;
};

typedef BackingStore<FieldStarLikelihoodRecord> FieldStarLikelihoodBackingStore;


struct SampleMassRecord
{
    Iteration iter;
    std::string starId;

    double primaryMass;
    double massRatio;

    double clusterMembership;
};

typedef BackingStore<std::vector<SampleMassRecord>> SampleMassBackingStore;

struct SampleWdMassRecord
{
    Iteration iter;
    std::string starId;

    double mass;
    double clusterMembership;
    double precursorLogAge;
    double coolingAge;
    double logTeff;
    double logLittleG;
};

typedef BackingStore<std::vector<SampleWdMassRecord>> SampleWdMassBackingStore;

#endif
