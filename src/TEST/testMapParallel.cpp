// This is simply going to test using maps in parallel
// Each worker will get a vector of ints
// The worker will create a map which counts the number of occurrences of each int
// It will then return the number of ints that occurred just once.
#include <vector>
#include <map>
#include <cstdlib>
#include <sstream>

#include "ProcessFramework2.h"

using namespace std;
typedef vector<int> IntVec;

IntVec makeRandomVec(size_t n, int maxVal)
{
    IntVec vec(n);
    for (size_t i = 0; i < n; i++)
        vec[i] = rand() % maxVal;
    return vec;
}

struct WorkItem
{
    WorkItem() {}
    WorkItem(const IntVec& vec) : vec_(vec) { }
    WorkItem(const WorkItem& wi) : vec_(wi.vec_) { }
    IntVec vec_;
};

struct Result
{
    Result(int val) : val_(val) {};
    int val_;
};

// Keeps generating copies of the same work item
class Generator
{
    public:
    Generator(size_t numItems, IntVec& vec) : numItems_(numItems), count_(0), vec_(vec) { };

    bool generate(WorkItem& item)
    {
        if (count_ >= numItems_) return false;
        item = WorkItem(vec_);
        count_++;
        return true;
    }

    size_t getNumConsumed() { return count_; }

    size_t numItems_; // number of work items to generate
    size_t count_;
    IntVec vec_;
};

// Does nothing
class PostProcessor
{
    public:
    PostProcessor() : count_(0)  { };
    void process(WorkItem& input, Result& result) { count_++; };
    size_t getNumProcessed() {return count_; };

    private:

    size_t count_;
};

class MapProcessor
{
    public:
    MapProcessor() { };
    Result process(WorkItem& item)
    {
        typedef map<int, int> CountMap;
        CountMap countMap;

        // create the count map
        for(IntVec::const_iterator i = item.vec_.begin();
            i != item.vec_.end();
            i++)
        {
            CountMap::iterator iter = countMap.find(*i);
            if (iter != countMap.end())
                iter->second++;
            else
                countMap.insert(CountMap::value_type(*i, 1));
        }

        // count the number of items with value 1
        Result res(0);
        for(CountMap::const_iterator iter = countMap.begin();
            iter != countMap.end();
            iter++)
        {
            if (iter->second == 1)
                res.val_++;
        }
        return res;
    }
};


void runNumThreads(int numThreads)
{
    size_t size = 1000;
    int maxVal = 500;
    size_t numItems = 10000;
    size_t bufferSize = 100;

    // Generate a random vector
    IntVec myVec = makeRandomVec(size, maxVal);

    // Create generator
    Generator myGenerator(numItems, myVec);

    // Create processors
    typedef MapProcessor Processor;
    vector<Processor *> processVec(numThreads);
    for(size_t i = 0; i < (size_t) numThreads; i++)
    {
        processVec[i] = new Processor();
    }

    // Create postprocessor
    PostProcessor postProcessor;
   
    // Run in parallel
    typedef ThreadScheduler<WorkItem, Result, Generator, Processor, PostProcessor> Scheduler;
    ostringstream oss;
    oss << "mapProcessor_" << numThreads;
    Scheduler scheduler(oss.str(), bufferSize);
    scheduler.processWorkParallel(myGenerator, processVec, &postProcessor);


    // Delete processors
    for(size_t i = 0; i < numThreads; i++)
    {
        delete processVec[i];
    }
};

int main()
{
    runNumThreads(1);
    runNumThreads(2);
    runNumThreads(4);
    runNumThreads(8);
}

