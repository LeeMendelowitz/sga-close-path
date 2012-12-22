// Some code to test the ThreadBase
#ifndef THREADDUMMY_H
#define THREADDUMMY_H

#include <iostream>
#include "ThreadBase.h"
#include <semaphore.h>

class DummyData
{
    public:
    std::vector<int> input; // new input
    int * pVal; // put the current sum here
};

class ThreadDummy : public ThreadBase<DummyData>
{
    public:
    ThreadDummy(sem_t* pReadySem, int id) :
        ThreadBase<DummyData>(pReadySem),
        id_(id),
        sum_(0)
        {};

    ~ThreadDummy()
    {
        std::cout << "Thread " << id_ << " dead! sum: " << sum_ << std::endl;
    }

    private:
    // Get new inputs.
    // Share results.
    void exchange(DummyData& d)
    {
        input_.swap(d.input);
        *d.pVal =  sum_;
    }

    void doWork()
    {
        size_t repeats = 100;
        for (size_t j =0; j < repeats; j++)
        {
            for (int i=0; i<input_.size(); i++)
                sum_ += input_[i];
        }
    }

    int id_;
    int sum_;
    std::vector<int> input_;
};

#endif
