#ifndef PATHCOUNTER_H
#define PATHCOUNTER_H
//****************************************************************
// pathCounter.h
// Lee Mendelowitz <lmendelo@umiacs.umd.edu>

// Use dynamic programming to compute the number of paths between
// two vertex pairs

// Note: The hashers are implemented within namespace std::tr1.
// I'm not sure about the compatability of doing this.
//****************************************************************

#include "HashMap.h"
#include "SGUtil.h"
#include <functional>

////////////////////////////////////////////////////////////////////////
class PathMapKey;
namespace std
{
    namespace tr1{
        template<>
        class hash<PathMapKey>
        {
            public:
            size_t operator()(const PathMapKey& key) const;
        };
    };
};


////////////////////////////////////////////////////////////////////////
// hash_combine
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
    std::tr1::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

////////////////////////////////////////////////////////////////////////
// e1_ = SENSE means leaving 3' end of v1_
// e1_ = ANTISENSE means leaving 5' end of v1_
// e2_ = SENSE means entering 5' end of v2_
// e2_ = ANTISENSE means entering 3' end of v1_
class PathMapKey
{
    public:
    PathMapKey(const Vertex * v1, EdgeDir e1, const Vertex * v2, EdgeDir e2, int d) :
        v1_(v1), e1_(e1), v2_(v2), e2_(e2), d_(d) {};
    PathMapKey() : v1_(0), e1_(ED_SENSE), v2_(0), e2_(ED_SENSE), d_(0) {};
    const Vertex * v1_;
    EdgeDir e1_; // end leaving v1_. SENSE meas
    const Vertex * v2_;
    EdgeDir e2_; // end entering v2_;
    int d_;

    //friend std::hash<PathMapKey>::operator()(const PathMapKey& key);

    bool operator==(const PathMapKey& other) const
    {
        return ((v1_ == other.v1_) &&
               (e1_ == other.e1_) &&
               (v2_ == other.v2_) &&
               (e2_ == other.e2_) &&
               (d_ == other.d_));
    }

    // Create the key such that (v1,e1) in the key is less than (v2,e2)
    static PathMapKey makeKey(const Vertex* v1, EdgeDir e1, const Vertex * v2, EdgeDir e2, int d)
    {
        bool oneIsLess = True;
        if (v1 < v2)
            oneIsLess = True;
        else if (v1 > v2)
            oneIsLess = False;
        else if (e1 < e2)
            oneIsLess = True;
        else if (e2 > e1)
            oneIsLess = False;
        if (oneIsLess)
            return PathMapKey(v1, e1, v2, e2, d);
        else:
            return PathMapKey(v2, e2, v1, e1, d);
    }
};


////////////////////////////////////////////////////////////////////////
// Define the hash function
namespace std
{
    namespace tr1{
//    template <>
    size_t hash<PathMapKey>::operator()(const PathMapKey& key) const
    {
        size_t seed = 0;
        hash_combine(seed, key.v1_);
        hash_combine(seed, key.v2_);
        return seed;
    }
    };
};


////////////////////////////////////////////////////////////////////////
typedef std::tr1::hash<PathMapKey> PathMapKeyHasher;
typedef DenseHashMap<PathMapKey, size_t, PathMapKeyHasher> PathCountMap;
class PathFinder
{
    public:
    PathFinder(StringGraph* pGraph, int maxGap, int maxOL);

    void initCountMap();

    private:
    int maxGap_;
    int maxOL_;
    StringGraph* pGraph_;
    PathCountMap countMap_;
};

#endif
