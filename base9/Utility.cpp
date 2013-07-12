#include <functional>
#include <thread>
#include <vector>

namespace base
{
    namespace utility
    {
        void parallelFor(const unsigned int size, std::function<void(const unsigned int)> func)
        {
            const unsigned int cores     = std::thread::hardware_concurrency();
            const unsigned int nbThreads = size < cores ? size : cores;

            std::vector < std::thread > threads;

            for (unsigned int idThread = 0; idThread < nbThreads; idThread++) {
                auto threadFunc = [=, &threads]() {
                    for (unsigned int i=idThread; i<size; i+=nbThreads) {
                        func(i);
                    }
                };

                threads.push_back(std::thread(threadFunc));
            }
            for (auto &t : threads)
                t.join();
        }
    }
}
