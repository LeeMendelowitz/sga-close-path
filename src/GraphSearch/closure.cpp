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

inline bool twinEqual(const Edge* e1, const Edge* e2)
{
    return (e1 == e2->getTwin());
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
        (*i)->getTwin()->setColor(c);
    }
}

// Return true if this closure includes a decision node in its interior
// A decision node has degree > 1 on at least one end
bool Closure::isDecisionClosure() const
{
    if (this->m_edges.size() < 2)
        return false;

    // Skip the last edge in closure so we only explore the interior vertices
    const EdgePtrVec::const_iterator E = this->m_edges.end()-1;
    for(EdgePtrVec::const_iterator iter = m_edges.begin(); iter != E; ++iter)
    {
        const Vertex * v = (*iter)->getEnd();
        if ( (v->countEdges(ED_SENSE) > 1) || (v->countEdges(ED_ANTISENSE) > 1) )
        {
            return true;
        }
    }
    return false;
}

// Add vertices on the interior of the walk to pVec
void Closure::getInteriorVertices(VertexPtrVec * pVec) const
{
    if (this->m_edges.empty())
        return;
    const EdgePtrVec::const_iterator E = this->m_edges.end()-1;
    for(EdgePtrVec::const_iterator iter = m_edges.begin(); iter != E; ++iter)
        pVec->push_back((*iter)->getEnd());
}

bool Closure::sfxOverlap(size_t startInd, const Closure& other, bool otherIsReverse)
{

    assert(startInd < m_edges.size());
    size_t n = m_edges.size() - startInd;
    if (n > other.m_edges.size())
        return false;

    EdgePtrVec::const_iterator s = m_edges.begin() + startInd;
    EdgePtrVec::const_iterator e = m_edges.end();
    if (otherIsReverse)
    {
        //     <--------- other
        // |--------> this
        // Check if the reverse complement of other (i.e. reverse twin path) is equal to this
        return equal(s, e, other.m_edges.rbegin(), twinEqual);
    }
    else
    {
        //        |-------------> other
        // |-------------> this
        return equal(s, e, other.m_edges.begin());
    }
}


// endInd is exclusive, zero based
bool Closure::pfxOverlap(size_t endInd, const Closure& other, bool otherIsReverse)
{
    assert(endInd > 0);
    assert(endInd <= m_edges.size());

    size_t n = endInd;
    if (n > other.m_edges.size())
        return false;

    EdgePtrVec::const_iterator s = m_edges.begin();
    EdgePtrVec::const_iterator e = s + n;
    if (otherIsReverse)
    {
        // other  <--------|
        //           |-----------------> this
        return equal(s, e, other.m_edges.rend() - n, twinEqual);
    }
    else
    {
        // other |---------->
        // this      |----------->
        return equal(s, e, other.m_edges.end()-n);
    }
}


// Compute the sequence of the closure
// Only include the last d1max of the first vertex's sequence,
// and the first d2max of the last vertex's sequence.
string Closure::computeSeq() const
{
    string walkSeq = this->computeFullSeq();

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

// Compute the full sequence of the closure
string Closure::computeFullSeq() const
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
    return walkSeq;
}

// Store the closure for the ClosePathResult if it is unique
void ClosureDB::process(const ClosePathResult& res)
{
    if (res.walks.size() != 1)
        return;
    Closure * c = new Closure(res.bundle->id, res.walks[0], res.bundle->d1max, res.bundle->d2max);
    nonContained_.push_back(c);
}

void ClosureDB::writeDecisionClosures(std::ostream& os, bool writeContained) const
{
    ClosurePtrVec::const_iterator iter = nonContained_.begin();
    for(; iter != nonContained_.end(); iter++)
    {
        const Closure * c = *iter;
        if (!c->isDecisionClosure())
            continue;
        os << c->id_ << "\t"
           << "NumEdges:" << c->getNumEdges() << "\t"
           << "Gap:" << c->getEndToStartDistance() << "\t";
        c->printWithOL(os);
        os << "\n";
    }

    if (writeContained)
    {
        ClosurePtrVec::const_iterator iter = contained_.begin();
        for(; iter != contained_.end(); iter++)
        {
            const Closure * c = *iter;
            if (!c->isDecisionClosure())
                continue;
            os << c->id_ << "\t"
               << "NumEdges:" << c->getNumEdges() << "\t"
               << "Gap:" << c->getEndToStartDistance() << "\t";
            c->printWithOL(os);
            os << "\n";
        }
    }
}


bool closureVecCmp(const Closure* p1, const Closure* p2)
{
    return ((*p1) < (*p2));
}

void ClosureDB::filterContainments()
{
    // Sort the Closures in ascending order of the number of edges
    // in the walks
    sort(nonContained_.begin(), nonContained_.end(), closureVecCmp);
    contained_.clear();
    ClosurePtrVec nonContained;

    // Remove closures which are contained in others 
    ClosurePtrVec::const_iterator i = nonContained_.begin();
    const ClosurePtrVec::const_iterator e = nonContained_.end();
    for (; i != e; i++)
    {
        bool isContained = false;
        Closure * c1 = *i;
        for (ClosurePtrVec::const_iterator b = i + 1; b != e; b++)
        {
            if ((*b)->contains(*c1))
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
    const ClosurePtrVec::const_iterator E = nonContained_.end();
    for(ClosurePtrVec::const_iterator iter = nonContained_.begin();
        iter != E;
        iter++)
    {
        // This will modify the graph and color interior edges of a closure red
        addClosurePath(pGraph, **iter, &iv);
    } 
    }

    // Erase edges interior to a closure
    size_t numEdgesRemoved = pGraph->sweepEdges(GC_RED);
    pGraph->setColors(GC_WHITE);
    pGraph->validate();

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
    pGraph->validate();

    #if DEBUG > 0
    cout << "Added " << nonContained_.size() << " closure paths to the graph as new nodes.\n"
         << "Removed " << numEdgesRemoved << " interior closure edges.\n"
         << "Removed " << numVertRemoved << " interior vertices which became islands.\n";
    #endif
}

ClosureDB::~ClosureDB()
{

    // Destroy all contained and non-conained Closure instance
    for(ClosurePtrVec::iterator iter = contained_.begin();
        iter != contained_.end();
        iter++)
    {
        delete *iter;
    }
    contained_.clear();

    for(ClosurePtrVec::iterator iter = nonContained_.begin();
        iter != nonContained_.end();
        iter++)
    {
        delete *iter;
    }
    nonContained_.clear();

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


void ClosureDB::indexClosures()
{
    firstEdgeMap_.clear();
    lastEdgeMap_.clear();
    const ClosurePtrVec::const_iterator E = nonContained_.end();
    for (ClosurePtrVec::const_iterator i = nonContained_.begin();
         i != E;
         i++)
    {
        const Closure * c = *i;
        const Edge * firstEdge = c->getFirstEdge();
        const Edge * lastEdge = c->getLastEdge();
        firstEdgeMap_.insert(EdgeClosureMap::value_type(firstEdge, c));
        lastEdgeMap_.insert(EdgeClosureMap::value_type(lastEdge, c));
    }
}

void ClosureDB::findClosureOverlaps()
{
    cout << "Finding closure overlaps" << endl;
    indexClosures();
    typedef pair<EdgeClosureMap::iterator, EdgeClosureMap::iterator> IterRange;
    const ClosurePtrVec::const_iterator E = nonContained_.end();
    for (ClosurePtrVec::const_iterator i = nonContained_.begin();
         i != E;
         i++)
    {
        Closure * c = *i;

        const size_t numEdges = c->m_edges.size();
        IterRange ret;
        int olSize;

        for (size_t ind = 0; ind < numEdges; ind++)
        {

            const Edge* e = c->m_edges[ind];
            const Edge* e_twin = e->getTwin();

            // Case 1: suffix overlap, other forward
            // this: -------->
            // other:    --------->
            ret = firstEdgeMap_.equal_range(e);
            for(EdgeClosureMap::iterator iter = ret.first; iter != ret.second; iter++)
            {
                const Closure * other = iter->second;
                if ((c > other) || (c==other && ind==0)) continue;
                if(c->sfxOverlap(ind, *other, false))
                {
//                    cout << "Found sfx overlap! Other is forward.\n";
//                    cout << "this: "; c->printWithOL(cout); cout << "\n";
//                    cout << "other: "; other->printWithOL(cout); cout << "\n";
                }

                //overlaps_.push_back(Overlap(c, ind, numEdges, other, 0, numEdges - ind)
            }

            // Case 2: suffix overlap, other reverse
            // this: -------->
            // other:    <--------
            ret = lastEdgeMap_.equal_range(e_twin);
            for(EdgeClosureMap::iterator iter = ret.first; iter != ret.second; iter++)
            {
                const Closure * other = iter->second;
                if (c > other) continue;
                if(c->sfxOverlap(ind, *other, true))
                {
//                    cout << "Found sfx overlap! Other is reverse.\n"; 
//                    cout << "this: "; c->printWithOL(cout); cout << "\n";
//                    cout << "other: "; other->printWithOL(cout); cout << "\n";
                }
            }

            // Case 3: prefix overlap, other forward
            // this:      -------->
            // other: ------->
            ret = lastEdgeMap_.equal_range(e);
            for(EdgeClosureMap::iterator iter = ret.first; iter != ret.second; iter++)
            {
                const Closure * other = iter->second;
                if ((c > other) || (c==other && ind==numEdges-1)) continue;
                if(c->pfxOverlap(ind+1, *other, false))
                {
//                    cout << "Found pfx overlap! Other is forward.\n"; 
//                    cout << "this: "; c->printWithOL(cout); cout << "\n";
//                    cout << "other: "; other->printWithOL(cout); cout << "\n";
                }
            }

            // Case 4: prefix overlap, other reverse
            // this:      -------->
            // other: <-------
            ret = firstEdgeMap_.equal_range(e_twin);
            for(EdgeClosureMap::iterator iter = ret.first; iter != ret.second; iter++)
            {
                const Closure * other = iter->second;
                if (c > other) continue;
                if(c->pfxOverlap(ind+1, *other, true))
                {
//                    cout << "Found pfx overlap! Other is reverse.\n"; 
//                    cout << "this: "; c->printWithOL(cout); cout << "\n";
//                    cout << "other: "; other->printWithOL(cout); cout << "\n";
                }
            }
        }
    }
}


// Add closures to graph. Mark edges/nodes for removal (but do not remove).
void ClosureAlgorithms::addClosuresToGraph(StringGraph * pGraph, const ClosureVec& closures, bool removeInteriorNodes)
{

    VertexPtrVec interiorNodes;

    for (ClosureVec::const_iterator iter = closures.begin();
         iter != closures.end();
         iter++)
    {
        const Closure& closure = *iter;
        closure.getInteriorVertices(&interiorNodes);
        ClosureAlgorithms::addClosureToGraph(pGraph, closure);
    }

    // If interior nodes became islands, remove them.
    // If removeInteriorNodes is true, remove all interior nodes of the walk.
    for (VertexPtrVec::const_iterator iter = interiorNodes.begin();
         iter != interiorNodes.end();
         iter++)
    {
        if (removeInteriorNodes)
            (*iter)->setColor(GC_RED);
        else if ((*iter)->countEdges() == 0)
            (*iter)->setColor(GC_RED);
    }
}

void ClosureAlgorithms::addClosureToGraph(StringGraph * pGraph, const Closure& c)
{

    if (c.getNumEdges() == 0)
        return;

    // Create a vertex for the path
    string walkSeq = c.computeFullSeq();
    size_t walklen = walkSeq.size();
    string walkId = c.id_;

    // Add vertex and edges to graph
    Vertex* walkv = new(pGraph->getVertexAllocator()) Vertex(walkId, walkSeq);
    pGraph->addVertex(walkv);

    // Add the overlaps for the vertex to the graph
    Vertex* leftv = c.getStartVertex();
    Edge* lefte = c.getFirstEdge();
    Vertex* rightv = c.getLastVertex();
    Edge* righte = c.getLastEdge();
    
    // Sanity checks
    assert(lefte->getStart() == leftv);
    assert(righte->getEnd() == rightv);

    // Remodel edges at the prefix of the walk
    EdgePtrVec edgesToFix = leftv->getEdges(!lefte->getDir());
    for (EdgePtrVec::iterator iter = edgesToFix.begin();
    iter != edgesToFix.end();
    iter++)
    {
        Edge * e = *iter; // original edge out of leftv
        Vertex * v_other = e->getEnd();
        Edge * etwin = e->getTwin(); // original edge from v_other
        bool otherIsRC = (etwin->getDir() == ED_ANTISENSE);

        // Create the overlap with the prefix of the walk
        SeqCoord other_coord = etwin->getMatchCoord();
        SeqCoord orig_coord = e->getMatchCoord();
        int matchLen = orig_coord.length();
        SeqCoord walk_coord = SeqCoord(0, matchLen-1, walklen);
        assert(other_coord.length() == walk_coord.length());

        // Add the new edge
        Overlap new_ovl(walkId, walk_coord, v_other->getID(), other_coord, otherIsRC, 0);
        Edge * newEdge = SGAlgorithms::createEdgesFromOverlap(pGraph, new_ovl, false);
        assert(newEdge != NULL);

        // Mark the old edges for removal
        e->setColor(GC_RED);
        etwin->setColor(GC_RED);
    }

    // Remodel edges at the suffix of walk
    edgesToFix = rightv->getEdges(!righte->getDir());
    for (EdgePtrVec::iterator iter = edgesToFix.begin();
    iter != edgesToFix.end();
    iter++)
    {
        Edge * e = *iter; // original edge out of rightv
        Vertex * v_other = e->getEnd();
        Edge * etwin = e->getTwin(); // original edge from v_other
        bool otherIsRC = (etwin->getDir() == ED_SENSE);

        // Create the overlap with the prefix of the walk
        SeqCoord other_coord = etwin->getMatchCoord();
        SeqCoord orig_coord = e->getMatchCoord();
        int matchLen = orig_coord.length();
        SeqCoord walk_coord = SeqCoord(walklen-matchLen, walklen-1, walklen);
        assert(other_coord.length() == walk_coord.length());

        // Add the new edge
        Overlap new_ovl(walkId, walk_coord, v_other->getID(), other_coord, otherIsRC, 0);
        Edge * newEdge = SGAlgorithms::createEdgesFromOverlap(pGraph, new_ovl, false);
        assert(newEdge != NULL);

        // Mark the old edges for removal
        e->setColor(GC_RED);
        etwin->setColor(GC_RED);
    }

    // Mark edges interior to the closure for deletion
    c.colorInteriorEdges(GC_RED);

    // Mark the first and last vertex of the walk for removal
    leftv->setColor(GC_RED);
    rightv->setColor(GC_RED);
}
