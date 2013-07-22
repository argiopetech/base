#include <functional>
#include <thread>
#include <vector>

#include "Utility.hpp"

namespace base
{
    namespace utility
    {
        void ThreadPool::parallelFor(const unsigned int size, std::function<void(const unsigned int)> func)
        {
            int loopThreads = (size < nThreads) ? size : nThreads;

            for (unsigned int idThread = 0; idThread < loopThreads; idThread++) {
                auto threadFunc = [=]() {
                    for (unsigned int i=idThread; i<size; i+=nThreads) {
                        func(i);
                    }
                };

                threads.at(idThread).addJob(threadFunc);
            }
            for (auto &t : threads)
                t.join();
        }
    }
}
