#include <algorithm>
#include <iostream>
#include <string>

#include "closure.h"
#include "closePathProcess.h"
#include "SGAlgorithms.h"

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

void Closure::colorInteriorEdges(GraphColor c) const
{
    const EdgePtrVec::const_iterator E = this->m_edges.end();
    for (EdgePtrVec::const_iterator i = this->m_edges.begin();
         i != E;
         i++)
    {
        (*i)->setColor(c);
    }
}

// Add vertices on the interior of the walk to pVec
void Closure::getInteriorVertices(VertexPtrVec * pVec) const
{
    const EdgePtrVec::const_iterator E = this->m_edges.end()-1;
    for(EdgePtrVec::const_iterator iter = m_edges.begin(); iter != E; ++iter)
        pVec->push_back((*iter)->getEnd());
}


// Compute the sequence of the closure
// Only include the last d1max of the first vertex's sequence,
// and the first d2max of the last vertex's sequence.
string Closure::computeSeq() const
{
    string walkSeq = this->getString(SGWT_START_TO_END);

    // Note that the walkSeq always has the starting vertex's sequence
    // oriented forward. If the edge out of the first vertex is ED_ANTISENSE,
    // then the first vertex sequence comes at the end of the walkSeq.
    // We would like the walk sequence ordered precisely as given by the m_edges
    // vector.
    // Take the reverse complement of the walk sequence if necessary.
    if (this->m_edges[0]->getDir() == ED_ANTISENSE)
        walkSeq = reverseComplement(walkSeq);


    Vertex * pX = this->getStartVertex();
    Vertex * pY = this->getLastVertex();
    size_t lX = pX->getSeqLen();
    size_t lY = pY->getSeqLen();
    assert(this->d1max_ <= (int) lX);
    assert(this->d2max_ <= (int) lY);

    //Only include the last d1max_ bp of the first vertex
    size_t trimLeft = lX - this->d1max_;
    //Only include the first d2max bp of the last vertex
    size_t trimRight = lY - this->d2max_;

    int lClosure = walkSeq.size() - trimLeft - trimRight;
    assert(lClosure >= 0);
    string seqClosure = walkSeq.substr(trimLeft, lClosure);
    return seqClosure;
}

// Store the closure for the ClosePathResult if it is unique
void ClosureDB::process(const ClosePathResult& res)
{
    if (res.walks.size() != 1)
        return;
    Closure c(res.bundle->id, res.walks[0], res.bundle->d1max, res.bundle->d2max);
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
                /*
                #if DEBUG > 0
                cout << "*******************\n";
                cout << "CONTAINMENT: " << (*b).id_ << "\n";
                (*b).printWithOL(cout);
                cout << "\n";
                (*i).printWithOL(cout);
                cout << "\n";
                #endif
                */
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

void ClosureDB::addClosurePaths(StringGraph* pGraph)
{

    // Iterate over the closures and add these paths to the graph as new nodes
    // Remove the closure edges from the graph
    // If any interior vertices are now islands, remove them from the graph
    pGraph->setColors(GC_WHITE);
    VertexPtrVec iv; // interior vertices

    {
    const ClosureVec::const_iterator E = nonContained_.end();
    for(ClosureVec::const_iterator iter = nonContained_.begin();
        iter != E;
        iter++)
    {
        // This will modify the graph and color interior edges of a closure red
        addClosurePath(pGraph, *iter, &iv);
    } 
    }

    // Erase edges interior to a closure
    size_t numEdgesRemoved = pGraph->sweepEdges(GC_RED);
    pGraph->setColors(GC_WHITE);

    // Get the unique list of interior vertices
    sort(iv.begin(), iv.end());
    iv.erase(unique(iv.begin(), iv.end()), iv.end());

    // Remove vertices which are interior to a closure which are now islands.
    size_t numVertRemoved = 0;
    const VertexPtrVec::const_iterator E = iv.end();
    for(VertexPtrVec::const_iterator iter = iv.begin();
        iter != E;
        iter++)
    {
        Vertex * pVertex = *iter;
        if (pVertex->countEdges()==0)
        {
            pGraph->removeIslandVertex(pVertex);
            numVertRemoved++;
        }
    }

    #if DEBUG > 0
    cout << "Added " << nonContained_.size() << " closure paths to the graph as new nodes.\n"
         << "Removed " << numEdgesRemoved << " interior closure edges.\n"
         << "Removed " << numVertRemoved << " interior vertices which became islands.\n";
    #endif
}


// Add the sequence of the closure path to the graph as a new node,
// with a single edge on either side.
// Delete edges within the path.
void ClosureDB::addClosurePath(StringGraph* pGraph, const Closure& c, VertexPtrVec* pInteriorVertices)
{
    // Create a vertex for the path
    string walkSeq = c.computeSeq();
    string walkId = c.id_;

    // Add the overlaps for the vertex to the graph
    Vertex* leftv = c.getStartVertex();
    Edge* lefte = c.getFirstEdge();
    Vertex* rightv = c.getLastVertex();
    Edge* righte = c.getLastEdge();
    
    // Sanity checks
    assert(lefte->getStart() == leftv);
    assert(righte->getEnd() == rightv);

    bool leftIsRc = (lefte->getDir() == ED_ANTISENSE);
    size_t leftv_len = leftv->getSeqLen();

    bool rightIsRc = (righte->getTwin()->getDir() == ED_SENSE);
    size_t rightv_len = rightv->getSeqLen();

    // Construct left overlap
    //cout << "s: " << leftv_len - c.d1max_ << " e: " << leftv_len-1 << " l: " << leftv_len << endl;

    // To avoid containment of the left or right node, trim the walkSequence by one base on either side if necessary.
    int leftOL = c.d1max_;
    int rightOL = c.d2max_;
    if (leftOL == (int) leftv_len)
    {
        leftOL--;
        walkSeq = walkSeq.substr(1,walkSeq.size()-1);
    }
    if (rightOL == (int) rightv_len)
    {
        rightOL--;
        walkSeq = walkSeq.substr(0,walkSeq.size()-1);
    }

    // Construct the coordinates of the overlaps
    SeqCoord lc(leftv_len - leftOL, leftv_len-1, leftv_len);
    SeqCoord wlc(0, leftOL-1, walkSeq.size());
    SeqCoord rc(0, rightOL-1, rightv_len);
    SeqCoord wrc(walkSeq.size()-rightOL, walkSeq.size()-1, walkSeq.size());
    if (leftIsRc) lc.flip();
    if (rightIsRc) rc.flip();
    std::string leftOlSeq = lc.getSubstring(leftv->getStr());
    std::string rightOlSeq = rc.getSubstring(rightv->getStr());
    if (leftIsRc) leftOlSeq = reverseComplement(leftOlSeq);
    if (rightIsRc) rightOlSeq = reverseComplement(rightOlSeq);
    cout << "-------------------------------------\n";
    cout << "left:          " << lc << endl;
    cout << "walk_left:     " << wlc << endl;
    cout << "right:         " << rc << endl;
    cout << "walk_right:    " << wrc << endl;
    cout << "leftv ol seq:  " << leftOlSeq << endl;
    cout << "walkv ol seq:  " << wlc.getSubstring(walkSeq) << endl;
    cout << "rightv ol seq: " << rightOlSeq << endl;
    cout << "walkv ol seq:  " << wrc.getSubstring(walkSeq) << endl;
    assert(lc.length() == wlc.length());
    assert(rc.length() == wrc.length());
    assert(leftOlSeq == wlc.getSubstring(walkSeq));
    assert(rightOlSeq == wrc.getSubstring(walkSeq));

    // Construct overlaps
    Overlap oleft(walkId, wlc, leftv->getID(), lc, leftIsRc, 0);  
    Overlap oright(walkId, wrc, rightv->getID(), rc, rightIsRc, 0);

    // Add vertex and edges to graph
    Vertex* walkv = new(pGraph->getVertexAllocator()) Vertex(walkId, walkSeq);
    pGraph->addVertex(walkv);
    Edge * newEdge = SGAlgorithms::createEdgesFromOverlap(pGraph, oleft, false);
    assert(newEdge != NULL);
    newEdge = SGAlgorithms::createEdgesFromOverlap(pGraph, oright, false);
    assert(newEdge != NULL);

    // Mark edges interior to the closure for deletion
    c.colorInteriorEdges(GC_RED);

    // Add interior vertices to pInteriorVertices
    c.getInteriorVertices(pInteriorVertices);
}
