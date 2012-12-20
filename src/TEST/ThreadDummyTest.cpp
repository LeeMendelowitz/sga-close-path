//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//
// Modified by Lee Mendelowitz (lmendelo@umiacs.umd.edu)
//-----------------------------------------------
//
// ProcessFramework - Generic framework for performing
// some operations on input data produced by a generator,
// serially or in parallel. 

#include "ThreadDummy.h"

#include <semaphore.h>
#include <pthread.h>
#include <cassert>

void test(int numThreads, int maxNums)
{

    // Helpful typedefs
    typedef std::vector<sem_t*> SemaphorePtrVector;
    typedef std::vector<ThreadDummy *> ThreadPtrVector;
    ThreadPtrVector threadVec(numThreads);
    SemaphorePtrVector semVec(numThreads);
    // Create the threads

    sem_t * threadAvailableSem = new sem_t;
    int ret = sem_init(threadAvailableSem, PTHREAD_PROCESS_PRIVATE, 0 );
    if(ret != 0)
    {
        std::cerr << "Semaphore initialization failed with error " << ret << "\n";
        std::cerr << "You are probably running on OSX which does not provide unnamed semaphores\n";
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < numThreads; ++i)
    {
        // Create and start the thread
        threadVec[i] = new ThreadDummy(threadAvailableSem, i);
        threadVec[i]->start();

        // Wait for the thread to become available
    }

    std::vector<int> intVec;
    intVec.reserve(maxNums);
    for (int i=0; i < maxNums; i++)
        intVec.push_back(i);

    std::vector<int> curVal(numThreads, 0);


    typedef std::vector<DummyData> WorkVector;
    WorkVector workVec;
    for (int i=0; i < numThreads; i++)
    {
        DummyData data;
        data.input = intVec;
        data.pVal = &curVal[i];
        workVec.push_back(data);
    }

    const int numRounds = 100;
    int round = 0;
    while(round < numRounds)
    {
        round++;
        // Give each thread a new number
        std::cout << "Waiting for a thread to become available." << std::endl;
        sem_wait(threadAvailableSem);
        int numReady = 0;
        for(int i = 0; i < numThreads; ++i)
        {
            ThreadDummy* pThread = threadVec[i];
            if (!pThread->isReady()) continue;
            numReady++;
            std::cout << "giving work " << " to thread " << i << std::endl;
            DummyData& d = workVec[i];
            pThread->exchangeData(d);
            std::cout << "Thread " << i << " reports partial sum " << curVal[i] << std::endl;
        }
        assert(numReady >= 1);
        for(int i = 1; i < numReady; i++)
        {
            int ret = sem_trywait(threadAvailableSem);
            assert(ret == 0);
        }
    }

    // Cleanup
    std::cout << " Stopping all threads..." << std::endl;
    for(int i = 0; i < numThreads; ++i)
    {
        std::cout << "stopping thread " << i << std::endl;
        threadVec[i]->stop(); // Blocks until the thread joins
    }
    std::cout << " Deleting all threads..." << std::endl;
    for(int i = 0; i < numThreads; ++i)
    {
        std::cout << "deleting thread " << i << std::endl;
        delete threadVec[i];
    }
    sem_destroy(threadAvailableSem);
}

int main()
{
    test(5, 100);
    return 1;
}
