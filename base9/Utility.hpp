#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <cassert>
#include <condition_variable>
#include <fstream>
#include <functional>
#include <ostream>
#include <iomanip>
#include <mutex>
#include <queue>
#include <thread>


#include "Settings.hpp"
#include "Star.hpp"

namespace base
{
    namespace utility
    {
        class ThreadPool
        {
          private:
            unsigned int numThreads(unsigned int threads)
            {
                return (threads > std::thread::hardware_concurrency()) ? std::thread::hardware_concurrency() : threads;
            }

            class WorkerThread
            {
              private:
                std::queue<std::function<void()>> workQueue;
                std::mutex m;
                std::condition_variable cv;

                std::unique_ptr<std::thread> localThread;

                bool joinable = false;
                bool terminate = false;

              public:
                WorkerThread()
                    : workQueue()
                {
                    localThread = std::unique_ptr<std::thread>(new std::thread([this]()
                                      {
                                          for(;;)
                                          {
                                              std::unique_lock<std::mutex> lk(m);

                                              while (this->workQueue.empty() && !terminate)
                                              {
                                                  joinable = false; // Can't join if we aren't running
                                                  cv.notify_one(); // Tell join it can finish
                                                  this->cv.wait(lk); // Wait for a signal
                                              }

                                              if (!terminate && !this->workQueue.empty()) // If we're still running...
                                              {
                                                  auto func = this->workQueue.front();
                                                  lk.unlock(); // Don't need lock while we run the function

                                                  func(); // It's a function pointer, so we run it with ()

                                                  lk.lock(); // Re-lock so we can mess with the queue
                                                  this->workQueue.pop();
                                              }
                                              else
                                              {
                                                  break;
                                              }
                                          }
                                      }));
                }

                ~WorkerThread()
                {
                    std::unique_lock<std::mutex> lk(m);
                    terminate = true; // Set the terminate flag
                    joinable = false; // Can't join, we're shutting down
                    cv.notify_all(); // Let everyone know we're ready to leave
                    lk.unlock(); // So that cv.wait() can pick it back up.

                    localThread->join(); // Wait for the thread to terminate on its own
                }

                void addJob(std::function<void()> f)
                {
                    std::lock_guard<std::mutex> lk(m); // Lock everything (while in scope)
                    joinable = true; // We have a job, therefore we're joinable
                    workQueue.push(f); // Push the job the queue
                    cv.notify_one(); // Signal the worker thread
                }

                void join()
                {
                    std::unique_lock<std::mutex> lk(m);
                    while (joinable) // Falls through if we're not joinable, sticks around till the queue is empty
                    {
                        cv.wait(lk);
                    }
                }
            };

            const unsigned int nThreads;

            std::vector < WorkerThread > threads;

          public:
            ThreadPool()
                : ThreadPool(std::thread::hardware_concurrency())
            {}

            ThreadPool(unsigned int nThreads)
                : nThreads(numThreads(nThreads)), threads(numThreads(nThreads))
            {}

            ThreadPool(const ThreadPool&) = delete;
            ThreadPool(const ThreadPool&&) = delete;
            ThreadPool& operator=(const ThreadPool&) = delete;

            void parallelFor(const unsigned int size, std::function<void(const unsigned int)> func);
        };

        std::pair<std::vector<std::string>, std::vector<StellarSystem>> readPhotometry (std::ifstream &fin, std::vector<double> &filterPriorMin, std::vector<double> &filterPriorMax, const Settings &settings);

        class InvalidPhotometry : public std::runtime_error
        {
            // using std::domain_error::domain_error;

          public:
            explicit InvalidPhotometry (const std::string& what_arg)
                : std::runtime_error(what_arg) {}

            explicit InvalidPhotometry (const char* what_arg)
                : std::runtime_error(what_arg) {}
        };

        std::ostream& format (std::ostream&);
    }
}

#endif
