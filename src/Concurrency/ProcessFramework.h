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
#include "ThreadWorker.h"
#include "Timer.h"

#ifndef PROCESSFRAMEWORK_H
#define PROCESSFRAMEWORK_H

template <class Input, class Output, class Generator, class Processor, class PostProcessor>
class ProcessFramework
{

    public:

    ProcessFramework(const std::string& name, size_t bufferSize = 1000, const size_t reportInterval = 1000) :
        name_(name),
        bufferSize_(bufferSize),
        reportInterval_(reportInterval)
    { };

    // Generic function to process n work items from a file. 
    // With the default value of -1, n becomes the largest value representable for
    // a size_t and all values will be read
//    template<class Input, class Output, class Generator, class Processor, class PostProcessor>
    size_t processWorkSerial(Generator& generator, Processor* pProcessor, PostProcessor* pPostProcessor, size_t n = -1)
    {
        Timer timer(name_, true);
        Input workItem;
        
        // Generate work items using the generic generation class while the number
        // of items consumed from the generator is less than n and there 
        // are still items to consume from the generator

        while(generator.getNumConsumed() < n && generator.generate(workItem))
        {
            Output output = pProcessor->process(workItem);
            
            pPostProcessor->process(workItem, output);
            if(generator.getNumConsumed() % reportInterval_ == 0)
                printf("[sga %s] Processed %zu items (%lfs elapsed)\n", name_.c_str(), generator.getNumConsumed(), timer.getElapsedWallTime());
        }

        assert(n == (size_t)-1 || generator.getNumConsumed() == n);

        //
        double proc_time_secs = timer.getElapsedWallTime();
        printf("[sga %s] processed %zu items in %lfs (%lf items/s)\n", 
                name_.c_str(), generator.getNumConsumed(), proc_time_secs, (double)generator.getNumConsumed() / proc_time_secs);    
        
        return generator.getNumConsumed();
    }

    // Design:
    // This function is a generic function to read some INPUT from a 
    // generic generator object, then perform work on them.
    // The actual processing is done by the Processor class 
    // that is passed in. The number of threads
    // created is determined by the size of the vector of processors - 
    // one thread per processor. 
    //
    // The function buffers batches of input data.
    // Once the buffers are full, the reads are dispatched to the thread
    // which run the actual processing independently. An optional post processor
    // can be specified to process the results that the threads return. If the n
    // parameter is used, at most n items will be read from the file
 //   template<class Input, class Output, class Generator, class Processor, class PostProcessor>
    size_t processWorkParallel(Generator& generator, 
                               std::vector<Processor*> processPtrVector, 
                               PostProcessor* pPostProcessor, 
                               size_t n = -1)
    {
        Timer timer(name_, true);

        // Helpful typedefs
        typedef ThreadWorker<Input, Output, Processor> Thread;
        typedef std::vector<Thread*> ThreadPtrVector;

        typedef std::vector<Input> InputItemVector;
        typedef std::vector<InputItemVector*> InputBufferVector;

        typedef std::vector<Output> OutputVector;
        typedef std::vector<OutputVector*> OutputBufferVector;
        typedef std::vector<sem_t*> SemaphorePtrVector;


        // Initialize threads, one thread per processor that was passed in
        int numThreads = processPtrVector.size();

        ThreadPtrVector threadVec(numThreads);
        InputBufferVector inputBuffers(numThreads);
        OutputBufferVector outputBuffers(numThreads);
        SemaphorePtrVector semVec(numThreads);

        // Create the threads
        for(int i = 0; i < numThreads; ++i)
        {
            semVec[i] = new sem_t;
            int ret = sem_init( semVec[i], PTHREAD_PROCESS_PRIVATE, 0 );
            if(ret != 0)
            {
                std::cerr << "Semaphore initialization failed with error " << ret << "\n";
                std::cerr << "You are probably running on OSX which does not provide unnamed semaphores\n";
                exit(EXIT_FAILURE);
            }

            // Create and start the thread
            threadVec[i] = new Thread(semVec[i], processPtrVector[i], bufferSize_);
            threadVec[i]->start();

            inputBuffers[i] = new InputItemVector;
            inputBuffers[i]->reserve(bufferSize_);

            outputBuffers[i] = new OutputVector;
            outputBuffers[i]->reserve(bufferSize_);
        }

        size_t numWorkItemsRead = 0;
        size_t numWorkItemsWrote = 0;
        bool done = false;
        int next_thread = 0;
        int num_buffers_full = 0;

        size_t reportInterval =  bufferSize_ * numThreads;
        if (reportInterval < reportInterval_)
        {
            reportInterval = reportInterval * (reportInterval_/reportInterval); 
        }

        while(!done)
        {
            // Parse reads from the stream and add them into the incoming buffers
            Input workItem;
            bool valid = generator.generate(workItem);
            if(valid)
            {
                inputBuffers[next_thread]->push_back(workItem);
                numWorkItemsRead += 1;

                // Change buffers if this one is full
                if(inputBuffers[next_thread]->size() == bufferSize_)
                {
                    ++num_buffers_full;
                    ++next_thread;
                }
            }
            
            done = !valid || generator.getNumConsumed() == n;

            // Once all buffers are full or the input is finished, dispatch the reads to the threads
            // by swapping work buffers. 
            if(num_buffers_full == numThreads || done)
            {
                int numLoops = 0;
                do
                {
                    // Wait for all threads to be ready to receive
                    for(int i = 0; i < numThreads; ++i)
                    {
                        sem_wait(semVec[i]);
                        Thread* pThread = threadVec[i];
                        pThread->swapBuffers(*inputBuffers[i], *outputBuffers[i]);
                    }
                    num_buffers_full = 0;
                    next_thread = 0;

                    // Process the results and clear the buffers
                    for(int i = 0; i < numThreads; ++i)
                    {
                        assert(inputBuffers[i]->size() == outputBuffers[i]->size());
                        for(size_t j = 0; j < inputBuffers[i]->size(); ++j)
                        {
                            pPostProcessor->process((*inputBuffers[i])[j], (*outputBuffers[i])[j]);
                            ++numWorkItemsWrote;
                        }
                        
                        inputBuffers[i]->clear();
                        outputBuffers[i]->clear();
                    }

                    if(generator.getNumConsumed() % (reportInterval) == 0)
                        printf("[sga %s] Processed %zu items\n", name_.c_str(), generator.getNumConsumed());

                    // This should never loop more than twice
                    assert(numLoops < 2);
                    ++numLoops;
                } while(done && numWorkItemsWrote < numWorkItemsRead);
            }
        }

        // Cleanup
        for(int i = 0; i < numThreads; ++i)
        {
            threadVec[i]->stop(); // Blocks until the thread joins
            delete threadVec[i];

            sem_destroy(semVec[i]);
            delete semVec[i];

            assert(inputBuffers[i]->empty());
            delete inputBuffers[i];

            assert(outputBuffers[i]->empty());
            delete outputBuffers[i];
        }
        assert(n == (size_t)-1 || generator.getNumConsumed() == n);
        assert(numWorkItemsRead == numWorkItemsWrote);

        double proc_time_secs = timer.getElapsedWallTime();
        printf("[sga %s] processed %zu items in %lfs (%lf items/s)\n", 
                name_.c_str(), generator.getNumConsumed(), proc_time_secs, (double)generator.getNumConsumed() / proc_time_secs);
        return generator.getNumConsumed();
    }


    private:
        const std::string name_;
        const size_t bufferSize_;
        const size_t reportInterval_;

};

#endif
