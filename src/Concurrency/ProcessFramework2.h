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
#ifndef PROCESSFRAMEWORK2_H
#define PROCESSFRAMEWORK2_H

//#define PROCESS_DEBUG 0

#include "ThreadBase.h"
#include "Timer.h"
#include <iostream>
#include <deque>

sem_t * makeSemaphore()
{
    sem_t * sem = NULL;
    sem = new sem_t;
    int ret = sem_init(sem, PTHREAD_PROCESS_PRIVATE, 0 );
    if(ret != 0)
    {
        std::cerr << "Semaphore initialization failed with error " << ret << "\n";
        std::cerr << "You are probably running on OSX which does not provide unnamed semaphores\n";
        exit(EXIT_FAILURE);
    }
    return sem;
}

#if PROCESS_DEBUG > 0
template <class T>
void printThreadStates(const std::vector<T*>& threadVec)
{
    for(size_t i = 0; i < threadVec.size(); i++)
    {
        const T * pThread = threadVec[i];
        std::cout << "Thread " << i << ": addr=" << pThread << " ready=" << (pThread ? pThread->isReady() : 0) << std::endl;
    }
}
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////
// Data type exchanged with a ProcessorThread
template <class Input, class Output>
struct ProcessorData
{
    typedef std::vector<Input> InputItemVector;
    typedef std::vector<Output> OutputItemVector;
    InputItemVector  inputItems;
    OutputItemVector outputItems;
};

// Data type exchanged with a GeneratorThread
template <class Input>
struct GeneratorData
{
    typedef std::vector<Input> InputItemVector;
    typedef std::vector<InputItemVector> InputBufferVector;
    InputBufferVector bufferVector;
};

// Data type exchanged with a PostProcessorThread
template <class Input, class Output>
struct PostProcessorData
{
    typedef ProcessorData<Input, Output> ProcData;
    typedef std::vector<ProcData> ProcDataVector;
    ProcDataVector items;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Define a Generator thread which fills a buffer of InputItemVectors, which will
// eventually get passed to a worker Thread.
template <class Input, class Generator>
class GeneratorThread : public ThreadBase< GeneratorData<Input> >
{

    typedef std::vector<Input> InputItemVector;
    typedef std::vector<InputItemVector> InputBufferVector;

    public:
    GeneratorThread(sem_t * pReadySemShared, Generator * pGenerator, size_t bufferSize, size_t numBuffers) :
        ThreadBase< GeneratorData<Input> >(pReadySemShared),
        pGenerator_(pGenerator),
        bufferSize_(bufferSize),
        numBuffers_(numBuffers),
        numGenerated_(0),
        generatorDone_(false),
        done_(false)

    { };

    bool isDone() const { return done_; }

    private:
    void exchange(GeneratorData<Input>& data);
    void setup();
    void doWork();

    Generator * pGenerator_;
    const size_t bufferSize_;
    const size_t numBuffers_;
    size_t numGenerated_;
    InputBufferVector inputBuffer_;
    volatile bool generatorDone_;
    volatile bool done_;

};


// Do a round of input generation at startup
template <class Input, class Generator>
void GeneratorThread<Input, Generator>::setup()
{
    doWork();
}

// Fill up the inputBuffer_ until it is of size numBuffers
template <class Input, class Generator>
void GeneratorThread<Input, Generator>::doWork()
{
    if (generatorDone_) return;

    assert(inputBuffer_.size() == 0);

    Input workItem;
    inputBuffer_.reserve(numBuffers_);
    for (size_t i = 0; i < numBuffers_ && !generatorDone_; i++)
    {
        // Generate items and place them into the items vector,
        // until either the generator has run out of items or the
        // items vector is full.
        inputBuffer_.push_back(InputItemVector());
        InputItemVector& items = inputBuffer_.back();
        items.reserve(bufferSize_);
        for (size_t j =0; j < bufferSize_; j++)
        {
            bool success = pGenerator_->generate(workItem);
            if (!success)
            {
                generatorDone_ = true;
                break;
            }
            items.push_back(workItem);
            numGenerated_++;
        }
    }
}

template <class Input, class Generator>
void GeneratorThread<Input, Generator>::exchange(GeneratorData<Input>& data)
{
    data.bufferVector.swap(inputBuffer_);
    inputBuffer_.clear();
    if (generatorDone_) done_ = true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Define a Post processor thread which consumes paired buffers of InputItemVectors and OutputItemVectors
template <class Input, class Output, class PostProcessor>
class PostProcessorThread : public ThreadBase< PostProcessorData<Input, Output> >
{

    typedef std::vector<Input> InputItemVector;
    typedef std::vector<Output> OutputItemVector;
    typedef std::vector<OutputItemVector> OutputBufferVector;
    typedef ProcessorData<Input, Output> ProcData;
    typedef std::vector<ProcData> ProcDataVector;

    public:
    PostProcessorThread(sem_t * pReadySemShared, PostProcessor * pPostProcessor) :
        ThreadBase< PostProcessorData<Input, Output> >(pReadySemShared),
        pPostProcessor_(pPostProcessor),
        numConsumed_(0)
    { };


    private:
    void exchange(PostProcessorData<Input, Output>& data);
    void doWork();

    PostProcessor * pPostProcessor_;
    ProcDataVector items_;
    size_t numConsumed_;

};

// Run the PostProcessor on all of the data
template <class Input, class Output, class PostProcessor>
void PostProcessorThread<Input, Output, PostProcessor>::doWork()
{
    const size_t numBuffers = items_.size();
    for (size_t i = 0; i < numBuffers; i++)
    {
        ProcData& procData = items_[i];
        InputItemVector& inputItems = procData.inputItems;
        OutputItemVector& outputItems = procData.outputItems;
        assert(inputItems.size() == outputItems.size());
        size_t numItems = inputItems.size();
        for (size_t j = 0; j < numItems; j++)
        {
            pPostProcessor_->process(inputItems[j], outputItems[j]);
            numConsumed_++;
        }
        inputItems.clear();
        outputItems.clear();
    }
    items_.clear();
}

template <class Input, class Output, class PostProcessor>
void PostProcessorThread<Input, Output, PostProcessor>::exchange(PostProcessorData<Input, Output>& data)
{

    // All work should be complete before accepting more
    assert(items_.empty());
    const size_t numItems = data.items.size();
    items_.reserve(numItems);
    for (size_t i = 0; i < numItems; i++)
    {
        items_.push_back(ProcData());
        ProcData& procData = items_.back();
        procData.inputItems.swap(data.items[i].inputItems);
        procData.outputItems.swap(data.items[i].outputItems);
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Define a Processor thread which performs work on an InputItemVector and populated an outputItemVector
template <class Input, class Output, class Processor>
class ProcessorThread : public ThreadBase< ProcessorData<Input, Output> >
{

    typedef std::vector<Input> InputItemVector;
    typedef std::vector<Output> OutputItemVector;
    typedef std::vector<OutputItemVector> OutputBufferVector;

    public:
    ProcessorThread(sem_t * pReadySemShared, Processor * pProcessor) :
        ThreadBase< ProcessorData<Input, Output> >(pReadySemShared),
        pProcessor_(pProcessor),
        numWorked_(0)
    { };


    private:
    void exchange(ProcessorData<Input, Output>& data);
    void doWork();

    Processor * pProcessor_;
    ProcessorData<Input, Output> data_;
    size_t numWorked_;

};

// Run the Processor on all of the data
template <class Input, class Output, class Processor>
void ProcessorThread<Input, Output, Processor>::doWork()
{
    size_t numItems = data_.inputItems.size();
    data_.outputItems.clear();
    data_.outputItems.reserve(numItems);
    for (size_t j = 0; j < numItems; j++)
    {
        Output output = pProcessor_->process(data_.inputItems[j]);
        data_.outputItems.push_back(output);
        numWorked_++;
    }
}

template <class Input, class Output, class Processor>
void ProcessorThread<Input, Output, Processor>::exchange(ProcessorData<Input, Output>& data)
{
    data_.inputItems.swap(data.inputItems); // Get new items to work
    data_.outputItems.swap(data.outputItems); // Give results
    data_.outputItems.clear();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////



template <class Input, class Output, class Generator, class Processor, class PostProcessor>
class ThreadScheduler
{

    public:

    ThreadScheduler(const std::string& name, size_t bufferSize = 1000, const size_t reportInterval = 1000) :
        name_(name),
        bufferSize_(bufferSize),
        reportInterval_(reportInterval)
    { };

    size_t processWorkSerial(Generator& generator, Processor* pProcessor, PostProcessor* pPostProcessor)
    {
        Timer timer(name_, true);
        Input workItem;
        
        // Generate work items using the generic generation class while the number
        // of items consumed from the generator is less than n and there 
        // are still items to consume from the generator
        while(generator.generate(workItem))
        {
            Output output = pProcessor->process(workItem);
            
            pPostProcessor->process(workItem, output);
            if(generator.getNumConsumed() % reportInterval_ == 0)
                printf("[sga %s] Processed %zu items (%lfs elapsed)\n", name_.c_str(), generator.getNumConsumed(), timer.getElapsedWallTime());
        }

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
    size_t processWorkParallel(Generator& generator, 
                               std::vector<Processor*> processPtrVector, 
                               PostProcessor* pPostProcessor)
    {
        Timer timer(name_, true);

        size_t consumed0 = generator.getNumConsumed();
        size_t processed0 = pPostProcessor->getNumProcessed();

        // Helpful typedefs

        typedef ProcessorThread<Input, Output, Processor> Worker;
        typedef std::vector<Worker *> ThreadPtrVector;

        typedef std::vector<Input> InputItemVector;
        typedef std::vector<InputItemVector*> InputBufferVector;

        typedef std::vector<Output> OutputVector;
        typedef std::vector<OutputVector*> OutputBufferVector;
        typedef std::vector<sem_t*> SemaphorePtrVector;

        typedef GeneratorThread<Input, Generator> GenThread;
        typedef PostProcessorThread<Input, Output, PostProcessor> PostProcThread;
        typedef ProcessorThread<Input, Output, Processor> ProcThread;

        typedef ProcessorData<Input, Output> ProcData;
        typedef GeneratorData<Input> GenData;
        typedef PostProcessorData<Input, Output> PostProcData;
        typedef std::vector<ProcData> ProcDataVector;

        const int numThreads = processPtrVector.size();

        // Create the Generator thread.
        sem_t * genSem = makeSemaphore();
        GenThread * pGenThread = new GenThread(genSem, &generator, bufferSize_, numThreads);
        pGenThread->start();

        // Create the Post Processor thread
        sem_t * postSem = makeSemaphore();
        PostProcThread * pPostProcThread = new PostProcThread(postSem, pPostProcessor);
        pPostProcThread->start();
        
        // Initialize worker threads, one thread per processor that was passed in
        ThreadPtrVector threadVec(numThreads);
        sem_t * sharedWorkerSem = makeSemaphore();

        // Create and start the worker threads
        for(int i = 0; i < numThreads; ++i)
        {
            // Create and start the thread
            threadVec[i] = new ProcThread(sharedWorkerSem, processPtrVector[i]);
            threadVec[i]->start();
        }

        // Set up data buffers
        std::deque<InputItemVector> inputDataQueue; // FIFO of input to be processed
        ProcDataVector processed; // Processed data buffer
        processed.reserve(2*numThreads);


        size_t numWorkItemsRead = 0;
        size_t numWorkItemsWrote = 0;
        size_t numWorkItemsProcessed = 0;
        size_t reportInterval = (reportInterval_/bufferSize_)*bufferSize_;
        if (reportInterval == 0) reportInterval = bufferSize_;
        size_t nextReport = reportInterval;


        // While there is still work to be generated or to be shared with the 
        // worker threads:
        bool done = false;
        while(!done)
        {
            // Get data from Generator, if necessary
            if (inputDataQueue.empty() && !pGenThread->isDone())
            {
                #if PROCESS_DEBUG > 0
                std::cout << "Getting input data from generator" << std::endl;
                #endif
                GenData data;
                sem_wait(genSem);
                pGenThread->exchangeData(data);
                size_t N = data.bufferVector.size();
                for(size_t i = 0; i < N; i++)
                {
                    inputDataQueue.push_back(InputItemVector());
                    InputItemVector & items = inputDataQueue.back();
                    items.swap(data.bufferVector[i]);
                    numWorkItemsRead += items.size();
                }
            }

            // Wait until a thread is ready for work. Give the first
            // ready thread some work.
            #if PROCESS_DEBUG > 0
            printThreadStates(threadVec);
            std::cout << "Waiting for thread" << std::endl;
            #endif

            assert(!inputDataQueue.empty());


            sem_wait(sharedWorkerSem); // this decrements the semaphore

            #if PROCESS_DEBUG > 0
            printThreadStates(threadVec);
            #endif

            bool foundReadyThread = false;
            for (size_t i = 0; i < (size_t) numThreads; i++)
            {
                ProcThread * pWorker = threadVec[i];
                if (!pWorker->isReady()) continue;


                #if PROCESS_DEBUG > 0
                std::cout << "Sending " << data.inputItems.size() << " items to thread " << i << std::endl;
                #endif

                // Create data to exchange with the thread
                ProcData data;
                data.inputItems.swap(inputDataQueue.front());
                inputDataQueue.pop_front();

                // Exchange data with the worker.
                pWorker->exchangeData(data);
                assert(data.inputItems.size() == data.outputItems.size());

                #if PROCESS_DEBUG > 0
                std::cout << "Received " << data.outputItems.size() << " items from thread " << i << std::endl;
                #endif

                numWorkItemsProcessed += data.outputItems.size();

                // If we received output, store it in the processed buffer.
                // Note: There will be no output after the first data exchange.
                if (data.inputItems.size() > 0)
                {
                    processed.push_back(ProcData());
                    ProcData& dest = processed.back();
                    data.inputItems.swap(dest.inputItems);
                    data.outputItems.swap(dest.outputItems);
                }

                foundReadyThread = true;
                break;
            }
            assert(foundReadyThread);


            // If the processed buffer is full, share with the post processor
            if(processed.size() >= (size_t) numThreads)
            {
                #if PROCESS_DEBUG > 0
                std::cout << "Sending output to post processor" << std::endl;
                #endif
                sem_wait(postSem);
                PostProcData data;
                data.items.swap(processed);
                for(size_t i = 0; i < data.items.size(); i++)
                    numWorkItemsWrote += data.items[i].outputItems.size();
                pPostProcThread->exchangeData(data);
                processed.clear();
            }

            if(numWorkItemsProcessed >= nextReport)
            {
                nextReport += reportInterval;
                printf("[sga %s] Processed %zu items (%lfs elapsed)\n", name_.c_str(), numWorkItemsProcessed, timer.getElapsedWallTime());
            }

            // If all the data has been distributed to the threads, break.
            done = pGenThread->isDone() && inputDataQueue.empty();
        }

        #if PROCESS_DEBUG > 0
        std::cout << "All work items have been given to worker threads" << std::endl;
        printThreadStates(threadVec);
        #endif


        // At this point all of the data has been distributed to the worker threads.
        // Wait for all threads to become available.
        ThreadPtrVector doneThreads;
        for (size_t i = 0; i < (size_t) numThreads; i++)
        {
            sem_wait(sharedWorkerSem);
        }

        // All the threads have now processed the remaining data.
        // Collect the data and give it to the post processor.
        #if PROCESS_DEBUG > 0
        printThreadStates(threadVec);
        #endif

        for (size_t i = 0; i < (size_t) numThreads; i++)
        {
            // Get this thread's data
            assert(threadVec[i]->isReady());
            ProcData data;
            threadVec[i]->exchangeData(data);
            assert(data.inputItems.size() == data.outputItems.size());

            // If we received output, store it in the processed buffer.
            if (data.inputItems.size() > 0)
            {
                processed.push_back(ProcData());
                ProcData& dest = processed.back();
                data.inputItems.swap(dest.inputItems);
                data.outputItems.swap(dest.outputItems);
            }
        }

        sem_wait(postSem);
        PostProcData data;
        data.items.swap(processed);
        for(size_t i = 0; i < data.items.size(); i++)
            numWorkItemsWrote += data.items[i].outputItems.size();
        pPostProcThread->exchangeData(data);
        processed.clear();

        // Delete the worker thread
        for(int i = 0; i < numThreads; ++i)
        {
            threadVec[i]->stop(); // Blocks until the thread joins
            delete threadVec[i];
        }

        // Delete the generator thread
        pGenThread->stop();
        delete pGenThread;
        pGenThread = NULL;

        // Delete the Post Processor thread
        sem_wait(postSem);
        pPostProcThread->stop();
        delete pPostProcThread;
        pPostProcThread = NULL;

        sem_destroy(genSem);
        delete genSem;

        sem_destroy(postSem);
        delete postSem;

        sem_destroy(sharedWorkerSem);
        delete sharedWorkerSem;

        // Check that the item counts are consistent
        assert(numWorkItemsRead == numWorkItemsWrote);
        assert(numWorkItemsRead == generator.getNumConsumed() - consumed0);
        assert(numWorkItemsWrote == pPostProcessor->getNumProcessed() - processed0);

        double proc_time_secs = timer.getElapsedWallTime();
        printf("[sga %s] processed %zu items in %lfs (%lf items/s)\n", 
                name_.c_str(), numWorkItemsRead, proc_time_secs, (double) numWorkItemsRead / proc_time_secs);
        return numWorkItemsRead;
    }


    private:
        const std::string name_;
        const size_t bufferSize_;
        const size_t reportInterval_;

};

#endif
