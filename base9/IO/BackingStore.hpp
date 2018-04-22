#ifndef IO_BACKING_STORE_HPP
#define IO_BACKING_STORE_HPP

#include <array>
#include <fstream>
#include <string>

#include "Cluster.hpp"
#include "constants.hpp"


enum class AdaptiveMcmcStage { FixedBurnin, AdaptiveBurnin, AdaptiveMainRun, MainRun };

struct Iteration
{
    const int val;
};

template <typename T>
class BackingStore
{
  public:
    virtual ~BackingStore() = default;

    virtual void save(Iteration, T) = 0;
    
    Iteration nextIteration() { return {++iteration}; }
    
  private:
    int iteration = 0;
};


template <typename T, typename HeaderData>
class FileBackingStore : public BackingStore<T>
{
  public:
    FileBackingStore(const std::string);

    virtual ~FileBackingStore();
        
  protected:
    std::ofstream fout;

    virtual void header(HeaderData) = 0;
};

#include "FileBackingStore.ipp"


template <typename T>
class SqlBackingStore : public BackingStore<T>
{
  public:
    virtual ~SqlBackingStore() = default;
};

#endif
