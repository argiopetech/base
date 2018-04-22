#ifndef IO_RECORDS_HPP
#define IO_RECORDS_HPP

#include "BackingStore.hpp"

struct MultiPopMcmcRecord
{
    const AdaptiveMcmcStage stage;

    const double lambda;

    const Cluster clust0;
    const Cluster clust1;
    
    const double logPost;
};

typedef BackingStore<MultiPopMcmcRecord> MultiPopBackingStore;


struct SinglePopMcmcRecord
{
    const AdaptiveMcmcStage stage;
    const Cluster clust;
    const double logPost;
};

typedef BackingStore<SinglePopMcmcRecord> SinglePopBackingStore;


struct StarParamsRecord
{
    const double fsLike;
    const std::vector<double> starData;
};

typedef BackingStore<StarParamsRecord> StarParamsBackingStore;

#endif
