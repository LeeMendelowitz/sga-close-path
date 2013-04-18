#include <algorithm>
#include "simplify.h"
#include <iostream>

#define DEBUG 1

using namespace std;

// Reverse an EdgePtrVec (using the edge twins in reverse order)
EdgePtrVec reverse(const EdgePtrVec& vec)
{
    // Initialze return vector
    EdgePtrVec vecOut(vec.size(), NULL);
    EdgePtrVec::iterator t = vecOut.begin();

    EdgePtrVec::const_reverse_iterator E = vec.rend();
    for (EdgePtrVec::const_reverse_iterator i = vec.rbegin();
         i != E;
         i++, t++)
        *t = (*i)->getTwin();
    return vecOut;
}

// Return true if vec1 contains vec2 as a subsequence
template< class T >
bool contains(const vector<T>& vec1, const vector<T>& vec2)
{
    size_t n1 = vec1.size();
    size_t n2 = vec2.size();
    if (n1 < n2)
        return false;
    typedef typename vector<T>::const_iterator Iter;
    const Iter b = vec2.begin();
    const Iter e = vec2.end();
    const Iter l = vec1.end() - n2;
    for (Iter i = vec1.begin(); 
         i <= l;
         i++)
    {
        if (std::equal(b, e, i))
            return true;
    }
    return false;
}


bool Closure::contains(const Closure& other) const
{
    const EdgePtrVec& vec1 = this->m_edges;
    const EdgePtrVec& vec2 = other.m_edges;
    if (vec1.size() < vec2.size())
        return false;

    if (::contains(vec1, vec2))
        return true;

    // See if vec1 contains the reverse of vec2
    EdgePtrVec vec2r = reverse(vec2);
    if (::contains(vec1, vec2r))
        return true;

    return false;
}

// Store the closure for the ClosePathResult if it is unique
void ClosureDB::process(const ClosePathResult& res)
{
    if (res.walks.size() != 1)
        return;
    Closure c(res.walks[0], res.bundle->d1max, res.bundle->d2max);
    nonContained_.push_back(c);
}

void ClosureDB::filterContainments()
{
    // Sort the Closures in ascending order of the number of edges
    // in the walks
    sort(nonContained_.begin(), nonContained_.end());
    contained_.clear();
    ClosureVec nonContained;

    // Remove closures which are contained in others 
    ClosureVec::const_iterator i = nonContained_.begin();
    const ClosureVec::const_iterator e = nonContained_.end();
    for (; i != e; i++)
    {
        bool isContained = false;
        const Closure& c1 = *i;
        for (ClosureVec::const_iterator b = i + 1; b != e; b++)
        {
            if ((*b).contains(c1))
            {
                isContained= true;
                break;
            }
        }

        if (isContained)
        {
            contained_.push_back(c1);
        }
        else
        {
            nonContained.push_back(c1);
        }
    }
    nonContained_ = nonContained;

    #if DEBUG > 0
    cout << "Marked " << contained_.size() << " as contained.\n"
         << "Marked " << nonContained_.size() << " as non contained.\n"
         << endl;

    #endif
}
