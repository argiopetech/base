#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <functional>

namespace base
{
    namespace utility
    {
        void parallelFor(const unsigned int size, std::function<void(const unsigned int)> func);
    }
}

#endif
