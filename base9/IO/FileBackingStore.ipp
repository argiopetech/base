#include <fstream>
#include <string>

#include "BackingStore.hpp"

using std::string;


template <typename T, typename HeaderData>
FileBackingStore<T, HeaderData>::FileBackingStore(const string filename)
    : fout(filename)
{
    if(!fout)
    {
        throw std::runtime_error(filename + " was not available for writing.");
    }
}


template <typename T, typename HeaderData>
FileBackingStore<T, HeaderData>::~FileBackingStore()
{
    fout.close();
}
