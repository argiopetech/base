#ifndef IO_BACKING_STORE_HPP
#define IO_BACKING_STORE_HPP

#include <array>
#include <fstream>
#include <functional>
#include <string>

#include "Cluster.hpp"
#include "constants.hpp"


enum class AdaptiveMcmcStage { FixedBurnin, AdaptiveBurnin, AdaptiveMainRun, MainRun };

struct RunData
{
    const signed long run;
};


struct Iteration
{
    const int val;
};


// TODO: Make this follow NVI
// http://www.gotw.ca/publications/mill18.htm
template <typename T>
class BackingStore
{
  public:
    virtual ~BackingStore() = default;

    virtual void save(T) = 0;

    virtual Iteration nextIteration() { return {++iteration}; }

  private:
    int iteration = 0;
};


// Provides a no-output version of the BackingStore which should match any FileBackingStore.
template <typename T>
class NullBackingStore : public BackingStore<T>
{
    public:
      virtual ~NullBackingStore() = default;

    virtual void save(T) {;}
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

#endif
