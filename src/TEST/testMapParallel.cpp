// This is simply going to test using maps in parallel
// Each worker will get a vector of ints
// The worker will create a map which counts the number of occurrences of each int
// It will then return the number of ints that occurred just once.
#include <vector>
#include <map>
#include <cstdlib>
#include <sstream>
#include <algorithm>

#define PROCESS_DEBUG 0
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
    WorkItem(const IntVec& vec) : vec_(&vec) { }
    WorkItem(const WorkItem& wi) : vec_(wi.vec_) { }
    const IntVec * vec_;
};

struct Result
{
    Result() {};
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
        for(IntVec::const_iterator i = item.vec_->begin();
            i != item.vec_->end();
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

class VectorProcessor
{
    public:
    VectorProcessor() { };
    Result process(WorkItem& item)
    {

        size_t N = item.vec_->size();
        //vec_.reserve(N);
        vec_.clear();
        vec_.insert(vec_.end(), item.vec_->begin(), item.vec_->end());
        //for(size_t i = 0; i< N; i++)
            //vec_[i] = item.vec_->operator[](i);

        sort(vec_.begin(), vec_.end());
        Result res(0);
        for (size_t i = 1; i < N-1; i++)
        {
            int v1 = vec_[i];
            int v0 = vec_[i-1];
            int v2 = vec_[i+1];
            if ((v0 != v1) && (v1 != v2))
            {
                res.val_++;
            }
        }
        if (vec_[0] != vec_[1]) res.val_++;
        if (vec_[N-1] != vec_[N-2]) res.val_++;

        int squareSum = 0;
        for (size_t i = 0; i < N; i++)
        {
            squareSum += vec_[i]; 
        }
        res.val_ += squareSum;
        return res;
    }

    IntVec vec_;
};

class VectorProcessorAlloc
{
    public:
    VectorProcessorAlloc() { };
    void process(WorkItem& item, Result& result)
    {
        result = Result(0);

        // Create a new vector of significant size.
        // This should require some memory allocation
        size_t N = 1024*1024*10/sizeof(int);
        std::vector<int *> myVec;
        myVec.reserve(N);
        for (size_t i = 0; i < N; i++)
        {
            int * p = new int;
            *p = i;
            myVec.push_back(p);
            result.val_ += i;
        }

        for (size_t i = 0; i < N; i++)
        {
            delete myVec[i];
        }
        return;
    }
};

template< class Processor>
void runNumThreads1(int numThreads, string& pfx)
{
    //size_t size = 1000000;
    size_t size = 100000;
    int maxVal = 500;
    //size_t bufferSize = 2000*4;
    size_t bufferSize = 2;
    size_t numItems = bufferSize*16*10;
    //size_t reportInterval = 16*bufferSize;
    size_t reportInterval = 20;

    // Generate a random vector
    IntVec myVec = makeRandomVec(size, maxVal);

    // Create generator
    Generator myGenerator(numItems, myVec);

    // Create processors
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
    oss << pfx << numThreads;
    Scheduler scheduler(oss.str(), bufferSize, reportInterval);
    scheduler.processWorkParallel(myGenerator, processVec, &postProcessor);


    // Delete processors
    for(size_t i = 0; i < numThreads; i++)
    {
        delete processVec[i];
    }
};

template< class Processor>
void runNumThreads2(int numThreads, string& pfx)
{
    //size_t size = 1000000;
    size_t size = 100000;
    int maxVal = 500;
    //size_t bufferSize = 2000*4;
    size_t bufferSize = 200;
    size_t numItems = bufferSize*16*10;
    //size_t reportInterval = 16*bufferSize;
    size_t reportInterval = 200;

    // Generate a random vector
    IntVec myVec = makeRandomVec(size, maxVal);

    // Create generator
    Generator myGenerator(numItems, myVec);

    // Create processors
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
    oss << pfx << numThreads;
    Scheduler scheduler(oss.str(), bufferSize, reportInterval);
    scheduler.processWorkParallel(myGenerator, processVec, &postProcessor);


    // Delete processors
    for(size_t i = 0; i < numThreads; i++)
    {
        delete processVec[i];
    }
};

int main()
{
    string pfx = "vecAlloc";
    //runNumThreads1<VectorProcessor>(1, pfx);
    runNumThreads1<VectorProcessorAlloc>(1, pfx);
    runNumThreads1<VectorProcessorAlloc>(2, pfx);
    runNumThreads1<VectorProcessorAlloc>(4, pfx);
    runNumThreads1<VectorProcessorAlloc>(8, pfx);
    runNumThreads1<VectorProcessorAlloc>(16, pfx);

    /*
    string pfx = "map";
    runNumThreads2<MapProcessor>(1, pfx);
    runNumThreads2<MapProcessor>(2, pfx);
    runNumThreads2<MapProcessor>(4, pfx);
    runNumThreads2<MapProcessor>(8, pfx);
    runNumThreads2<MapProcessor>(16, pfx);
    */

    /*
    pfx = "map";
    runNumThreads2<MapProcessor>(1, pfx);
    runNumThreads2<MapProcessor>(2, pfx);
    runNumThreads2<MapProcessor>(4, pfx);
    */
}

