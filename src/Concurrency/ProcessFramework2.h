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

#ifndef PROCESSFRAMEWORK2_H
#define PROCESSFRAMEWORK2_H

sem_t * makeSemaphore()
{
    sem_t * sem = NULL;
    sem = new sem_t;
    int ret = sem_init( sem, PTHREAD_PROCESS_PRIVATE, 0 );
    if(ret != 0)
    {
        std::cerr << "Semaphore initialization failed with error " << ret << "\n";
        std::cerr << "You are probably running on OSX which does not provide unnamed semaphores\n";
        exit(EXIT_FAILURE);
    }
    return sem;
}

void semGetValueFailed()
{
    if (ret != 0)
    {
        std::cerr << "Semaphore get_value failed with error " << ret << "\n";
        exit(EXIT_FAILURE);
    }
}

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
}:

// Data type exchanged with a PostProcessorThread
template <class Input, class Output>
struct PostProcessorData
{
    typedef ProcessorData<Input, Output> PairedData;
    typedef std::vector<PairedData> PairedDataVector;
    PairedDataVector pairedDataVector;
}:

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
        bool done_(false)

    { 
        inputBuffer_ = InputBufferVector(numBuffers);
    };

    private:
    Generator * pGenerator_;
    const size_t bufferSize;
    const size_t numBuffers;
    size_t numGenerated;
    InputBufferVector inputBuffer_;
    bool done_;

    void doWork();
    void exchange(GeneratorData<Input>& data);
};

// Fill up the inputBuffer_ until it is of size numBuffers
template <class Input, class Generator>
void GeneratorThread<Input, Generator>::doWork()
{
    if (done_) return;

    assert(inputBuffer_.size() == 0);

    Input workItem;
    inputBuffer_.reserve(numBuffers);
    for (size_t i = 0; i < numBuffers && !done_; i++)
    {
        // Generate items and place them into the items vector,
        // until either the generator has run out of items or the
        // items vector is full.
        inputBuffer_.push_back(InputItemVector())
        InputItemVector& items = inputBuffer_.back();
        items.reserve(bufferSize);
        while(items.size() < bufferSize)
        {
            bool success = pGenerator_->generate(workItem);
            if (!success)
            {
                done_ = true;
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
    data.inputBufferVector.swap(inputBuffer_);
    inputBuffer_.clear();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Define a Post processor thread which consumes paired buffers of InputItemVectors and OutputItemVectors
template <class Input, class Output, class PostProcessor>
class PostProcessorThread : public ThreadBase< PostProcessorData<Input, Output> >
{

    typedef std::vector<Input> InputItemVector;
    typedef std::vector<Output> OutputItemVector;
    typedef std::vector<Output> OutputItemVector;
    typedef std::vector<OutputItemVector> OutputBufferVector;
    typedef ProcessorData<Input, Output> PairedData;
    typedef std::vector<PairedData> PairedDataVector;

    public:
    PostProcessorThread(sem_t * pReadySemShared, PostProcessor * pPostProcessor) :
        ThreadBase< PostProcessorData<Input, Output> >(pReadySemShared),
        pPostProcessor_(pPostProcessor),
        numConsumed_(0)
    { }

    private:
    PostProcess * pPostProcessor_;
    PairedDataVector pairedData_;
    size_t numConsumed_;

    void doWork();
    void exchange(PostProcessorData<Input, Output>& data);
};

// Run the PostProcessor on all of the data
template <class Input, class Output, class PostProcessor>
void PostProcessorThread<Input, Output, PostProcessor>::doWork()
{
    size_t numBuffers = pairedData_.size();
    for (size_t i = 0; i < numBuffers; i++)
    {
        PairedData& pairedData = pairedData_[i];
        InputItemVector& inputItems = pairedData.inputItems;
        OutputItemVector& outputItems = pairedData.outputItems;
        assert(inputItems.size() == outputItems.size());
        size_t numItems = itemsItems.size();
        for (size_t j = 0; j < numItems; j++)
        {
            pPostProcessor->process(inputItems[j], outputItems[j]);
            numConsumed_++;
        }
        inputItems.clear();
        outputItems.clear();
    }
    pairedData_.clear();
}

template <class Input, class Output, class PostProcessor>
void PostProcessorThread<Input, Output, PostProcessor>::exchange(PostProcessorData<Input, Output>& data)
{
    data.pairedDataVector.swap(pairedData_);
    data.pairedDataVector.clear();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Define a Processor thread which performs work on an InputItemVector and populated an outputItemVector
template <class Input, class Output, class Processor>
class ProcessorThread : public ThreadBase< ProcessorData<Input, Output> >
{

    typedef std::vector<Input> InputItemVector;
    typedef std::vector<Output> OutputItemVector;
    typedef std::vector<Output> OutputItemVector;
    typedef std::vector<OutputItemVector> OutputBufferVector;

    public:
    ProcessorThread(sem_t * pReadySemShared, Processor * pProcessor) :
        ThreadBase< ProcessorData<Input, Output> >(pReadySemShared),
        pProcessor_(pProcessor),
        numWorked_(0)
    { }

    private:
    PostProcess * pPostProcessor_;
    ProcessorData<Input, Output> data_;
    size_t numWorked_;

    void doWork();
    void exchange(ProcessorData<Input, Output>& data);
};

// Run the Processor on all of the data
template <class Input, class Output, class Processor>
void ProcessorThread<Input, Output, Processor>::doWork()
{
    size_t numItems = data.inputItems.size();
    data.outputItems_.clear();
    data.outputItems_.reserve(numItems);
    for (size_t j = 0; j < numItems; j++)
    {
        Output output = pProcessor->process(workItem);
        data.outputItems_.push_back(output);
        numWorked_++;
    }
}

template <class Input, class Output, class Processor>
void ProcessorThread<Input, Output, Processor>::exchange(ProcessorData<Input, Output>& data)
{
    data_.inputItems_.swap(data.inputItems_); // Get new items to work
    data_.outputItems_.swap(data.outputItems_); // Give results
    data_.outputItems_.clear();
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

    // Generic function to process n work items from a file. 
    // With the default value of -1, n becomes the largest value representable for
    // a size_t and all values will be read
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

        const int numThreads = processPtrVector.size();

        // Create the Generator thread. Get the generator going by exchanging with it.
        sem_t * genSem = makeSemaphore();
        GenThread * pGenThread = new GenThread(genSem, &generator, bufferSize_, numThreads);
        GenData dummy;
        pGenThread->exchangeData(dummy);

        // Create the Post Processor thread
        sem_t * postSem = makeSemaphore();
        PostPocThread * pPostProcThread = new PostProcThread(postSem, pPostProcessor);

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


        // TO DO:
        // - Write code for exchanging with generator
        // - Wait until thread is available. Get Results, store for PostProcessor.
        // - Hand to PostProcessor
        // FIFO queue for storing InputItemVector's generated by the Generator Thread
        std::deque<InputItemVector> inputDataQueue;



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

            // TO DO: KEEP THE inputBuffers. Only fill those buffers that are empty.
            // with new values. The buffers should be emptied when they are shipped off for post-processing.
            // Get rid of the inputBuffer object.



            // If we have room in the inputBuffer, accept a new batch of InputItemVectors from
            // the generator thread.
            if (inputBuffer.size() + numGeneratorItems < MAX_BUFFER_VECTOR)
            {
                InputBufferVector newItemVectors;
                sem_wait(genSem);
                pGenWorker->swapBuffers(dummyGeneratorInput, newItemVectors);
                // Put the new inputItemVectors into the inputBuffer
                inputBuffer.insert(inputBuffer.end(), newItemVectors.begin(), newItemVectors.end());
            }

            done = !valid || generator.getNumConsumed() == n;

            // Once all buffers are full or the input is finished, dispatch the reads to the threads
            // by swapping work buffers. 
            if(num_buffers_full == numThreads || done)
            {
                int numLoops = 0;
                do
                {
                    // Wait for at least one thread to be ready to receive
                    sem_wait(workerSem);
                    int numReady = 0;

                    for(int i = 0; i < numThreads; ++i)
                    {
                        Thread* pThread = threadVec[i];
                        if (!pThread->isReady()) continue;
                        numReady++;

                        // Pop inputBuffers from the front
                        InputItemVector * pInput = inputBuffer.front();
                        inputBuffer.pop_front();
                        inputBuffers[i] = pInput;
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

                    if(generator.getNumConsumed() % reportInterval_ == 0)
                        printf("[sga %s] Processed %zu items (%lfs elapsed)\n", name_.c_str(), generator.getNumConsumed(), timer.getElapsedWallTime());
                    //if(generator.getNumConsumed() % (reportInterval) == 0)
                     //   printf("[sga %s] Processed %zu items\n", name_.c_str(), generator.getNumConsumed());

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
