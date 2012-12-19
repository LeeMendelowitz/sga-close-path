///-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//
// Modified by Lee Mendelowitz (lmendelo@umiacs.umd.ed)
//-----------------------------------------------
//
// ThreadBase - Generic thread base class.
// The ThreadBase has a 

// to perform a batch of work. It gets a vector of
// Input items from the master thread, uses Processor to
// perform some operation on the data which returns a value 
// of type Output. A vector of output is swapped back to the master
// thread
//
#ifndef THREADBASE_H
#define THREADBASE_H

#include <semaphore.h>
#include "Util.h"

template <class Data>
class ThreadBase
{

    public:
        ThreadBase(sem_t* pReadySemShared);
        ~ThreadBase();

        // This is called externally.
        // Accept new inputs and communicate outputs,
        // through the Data object.
        // The real exchange is done in exchange(), provided by the derived class.
        void exchangeData(Data& D);

        // External control functions
        void start();
        void stop();
        bool isReady();

    protected:

        // Main work loop
        void run();

        // Derived class should provide this function to
        // do the work.
        virtual void doWork() = 0;

        // Communicate results and accept new inputs.
        // Derived class must provide this function.
        virtual void exchange(Data& D) = 0;

        // Thread entry point
        static void* startThread(void* obj);
        
        // Handles
        pthread_t m_thread;

        // External semaphore to post to
        // when the thread is ready to receive data
        sem_t* m_pReadySemShared;

        // Internal semaphore to post to when
        // the thread is ready to receive data
        sem_t m_readySem;

        // Semaphore the external caller posts to when
        // data is ready to consume
        sem_t m_producedSem;
        
        // Shared data
        pthread_mutex_t m_mutex;

        volatile bool m_stopRequested;
        bool m_isReady;
};

// Implementation
template<class Data>
ThreadBase<Data>::ThreadBase(sem_t* pReadySemShared) :
                                                  m_pReadySemShared(pReadySemShared),
                                                  m_stopRequested(false), 
                                                  m_isReady(false)
{
    // Set up semaphores and mutexes
    int ret1 = sem_init(&m_readySem, PTHREAD_PROCESS_PRIVATE, 0 );
    int ret2 = sem_init(&m_producedSem, PTHREAD_PROCESS_PRIVATE, 0 );
    if( (ret1!=0) || (ret2!=0))
    {
        std::cerr << "Semaphore initialization failed with error " << ret1 << " and " << ret2 <<  "\n";
        std::cerr << "You are probably running on OSX which does not provide unnamed semaphores\n";
        exit(EXIT_FAILURE);
    }

    int ret = pthread_mutex_init(&m_mutex, NULL);
    if(ret != 0)
    {
        std::cerr << "Mutex initialization failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
}

//
template<class Data>
ThreadBase<Data>::~ThreadBase()
{
    sem_destroy(&m_producedSem);
    sem_destroy(&m_readySem);
    int ret = pthread_mutex_destroy(&m_mutex);
    if(ret != 0)
    {
        std::cerr << "Mutex destruction failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
    
}

// Externally-called function to start the worker
template<class Data>
void ThreadBase<Data>::start()
{
    int ret = pthread_create(&m_thread, 0, &ThreadBase<Data>::startThread, this);
    if(ret != 0)
    {
        std::cerr << "Thread creation failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
    
}

// Externally-called function to tell the worker to stop
template<class Data>
void ThreadBase<Data>::stop()
{
    // Wait for the work thread to finish
    // with the data its processing then set the stop
    // flag and unblock the work process
    sem_wait(&m_readySem);
    pthread_mutex_lock(&m_mutex);
    m_stopRequested = true;
    pthread_mutex_unlock(&m_mutex);
    sem_post(&m_producedSem);
    int ret = pthread_join(m_thread, NULL);
    if(ret != 0)
    {
        std::cerr << "Thread join failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
}

//
template<class Data>
void ThreadBase<Data>::exchangeData(Data& d)
{
    sem_wait(&m_readySem);
    pthread_mutex_lock(&m_mutex);
    m_isReady = false;
    storeData(d);
    sem_post(&m_producedSem);
    pthread_mutex_unlock(&m_mutex);
}

//
template<class Data>
void ThreadBase<Data>::run()
{
    // Indicate that the thread is ready to receive data
    pthread_mutex_lock(&m_mutex);
    m_isReady = true;
    sem_post(&m_readySem);
    sem_post(m_pReadySemShared);
    pthread_mutex_unlock(&m_mutex);

    while(1)
    {
        // Block until there is some data to use
        sem_wait(&m_producedSem);

        // If the stop flag is now set, finish
        if(m_stopRequested)
            break;

        pthread_mutex_lock(&m_mutex);

        // Do the work
        doWork();

        m_isReady = true;
        pthread_mutex_unlock(&m_mutex);

        // Post to the semaphore
        sem_post(m_pReadySemShared);
        sem_post(&m_readySem);
    }
}

// 
template<class Data>
bool ThreadBase<Data>::isReady()
{
    return m_isReady;
}


// Thread entry point
template<class Data>
void* ThreadBase<Data>::startThread(void* obj)
{
    reinterpret_cast<ThreadBase*>(obj)->run();
    return NULL;
}



#endif
