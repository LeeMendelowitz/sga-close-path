#include "contigMerger.h"
#include "closePathProcess.h"
#include "SGUtil.h"
#include "ScaffoldGraph.h"
#include "ScaffoldVisitors.h"
#include "ScaffoldAlgorithms.h"
#include "ScaffoldLink.h"
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std;

// For storing ClosureDatas in map
ClosureDataKey::ClosureDataKey(const ClosureData& c) :
        v1_(c.b_.vertex1ID),
        v2_(c.b_.vertex2ID),
        dir1_(c.b_.dir1),
        dir2_(c.b_.dir2)
{
    reorient();
}

ClosureDataKey::ClosureDataKey(const std::string& root, const ScaffoldLink& link)
{
    v1_ = root;
    v2_ = link.endpointID;
    dir1_ = link.getDir();
    dir2_ = link.getTwinDir();
    reorient();
}

void ClosureDataKey::reorient()
{
    if (v1_ > v2_)
    {
        swap(v1_, v2_);
        swap(dir1_, dir2_);
    }
}

std::ostream& operator<<(std::ostream& os, const ClosureDataKey& key)
{
    os << "v1: " << key.v1_ << " dir1: " << key.dir1_
       << " v2: " << key.v2_ << " dir2: " << key.dir2_;
    return os;
}

bool ClosureDataKey::operator<(const ClosureDataKey& other) const
{
    bool lt[4] = {v1_ < other.v1_, v2_ < other.v2_, dir1_ < other.dir1_, dir2_ < other.dir2_};
    bool gt[4] = {v1_ > other.v1_, v2_ > other.v2_, dir1_ > other.dir1_, dir2_ > other.dir2_};

    for (size_t i = 0; i < 4; i++)
    {
        if (lt[i])
            return true;
        else if (gt[i])
            return false;
    }
    return false;
}


// For sorting closures by the the vertexes they close
bool cmpClosureData(const ClosureData& c1, const ClosureData& c2)
{
    string c1v1,c1v2, c2v1, c2v2;
    c1v1 = c1.getVertex1();
    c1v2 = c1.getVertex2();
    if (c1v2 < c1v1) swap(c1v1, c1v2);
    c2v1 = c2.getVertex1();
    c2v2 = c2.getVertex2();
    if (c2v2 < c2v1) swap(c2v1, c2v2);

    if (c1v1 < c2v1)
        return true;
    else if (c1v1 > c2v1)
        return false;
    else if (c1v2 < c2v2)
        return true;
    else if (c1v2 > c2v2)
        return false;
    return true;
}

ContigMerger::ContigMerger(StringGraph* pGraph, const std::string& astatFile, float astatThreshold, int minLength) :
    pGraph_(pGraph),
    astatThreshold_(astatThreshold),
    minLength_(minLength)

{

    // Read astat records from the astat file
    AstatReader::readAstatFile(astatFile, astatMap_);

    // Populate the set of single copy contigs
    AstatMap::const_iterator iter = astatMap_.begin();
    const AstatMap::const_iterator E = astatMap_.end();
    int nb = 0;
    for (; iter != E; iter++)
    {
        const AstatRecord& rec = iter->second;
        if (rec.len >= minLength_ && rec.astat >= astatThreshold_)
        {
            singleCopy_.insert(rec.contig);
            nb += rec.len;
        }
    }
    
    cout << "Marked " << singleCopy_.size() << " contigs"
         << " (" << nb << " bp) as single copy.\n";
}

void ContigMerger::process(const ClosePathResult& res)
{

    const int MIN_LINKS = 2;

    // Only store results for closures between single copy contigs
    // Do not allow self closures
    if (singleCopy_.count(res.bundle->vertex1ID) == 1 &&
        singleCopy_.count(res.bundle->vertex2ID) == 1 &&
        res.walks.size() == 1 &&
        res.bundle->n >= MIN_LINKS)
    {
        ClosureData data(res.bundle->id, res.walks, *res.bundle);
        closures_.push_back(data);
        ClosureDataKey key(data);
        //cout << "Trying to add key: " << key << std::endl;
        assert(closureMap_.count(key)==0);
        closureMap_.insert(ClosureMap::value_type(key, data));
    }
}

void ContigMerger::postProcess()
{

    cout << "-------------------\n"
         << "CONTIG MERGER:\n";


    // Remove closures

    // summarize stored closures
    int numUnique = 0;
    for( ClosureVec::const_iterator iter = closures_.begin();
         iter != closures_.end();
         iter++)
    {
        if (iter->walks_.size() == 1)
            numUnique++;
        const Bundle& b = iter->b_;
        const std::string& id1 = b.vertex1ID;
        const std::string& id2 = b.vertex2ID;
        size_t l1 = pGraph_->getVertex(id1)->getSeqLen();
        size_t l2 = pGraph_->getVertex(id2)->getSeqLen();

        // Make a string of walk lengths
        ostringstream oss;
        oss << iter->walks_.size() << ":";
        for (SGWalkVector::const_iterator witer = iter->walks_.begin();
             witer != iter->walks_.end();
             witer++)
        {
            oss << "numEdges:" << witer->getNumEdges() << "_Dist:" << witer->getEndToStartDistance();
            if (witer != iter->walks_.end()-1)
                oss << ",";
        }
        string walkStr = oss.str();


        cout << setprecision(5);
        cout << id1 << "(" << b.dir1 << ") " << " L: " << l1 << " - "
             << id2 << "(" << b.dir2 << ") " << " L: " << l2
             << " n: " << b.n
             << " gapEst: " << b.gap << " "
             << " std: " << b.std << " "
             << walkStr << "\n";

    }

    cout << "Found " << closures_.size() << " bundles with closures between single copy contigs.\n";
    cout << "Found " << numUnique << " bundles with unique closures between single copy contigs.\n";
    cout << "DONE: CONTIG MERGER\n"
         << "--------------------\n";

    applyClosuresToGraph();
}

// Build a scaffold graph using closures between single copy contigs
void ContigMerger::applyClosuresToGraph()
{
    ScaffoldGraph scaffGraph;
    const int MIN_LINKS = 2;

    // Add vertices to the Scaffold Graph
    size_t numVerticedAdded = 0;
    size_t numEdgesAdded = 0;
    for(ClosureVec::const_iterator iter = closures_.begin();
        iter != closures_.end();
        iter++)
    {
        if (iter->walks_.size() != 1) continue;
        if (iter->b_.n < MIN_LINKS) continue;
        const Bundle& b = iter->b_;

        string v1 = iter->getVertex1();
        string v2 = iter->getVertex2();
        const Vertex* pV1 = pGraph_->getVertex(v1);
        const Vertex* pV2 = pGraph_->getVertex(v2);
        assert(pV1 != NULL);
        assert(pV2 != NULL);
        ScaffoldVertex * pScaffV1 = scaffGraph.getVertex(v1);
        ScaffoldVertex * pScaffV2 = scaffGraph.getVertex(v2);
        if (pScaffV1 == NULL)
        {
            pScaffV1 = new ScaffoldVertex(v1, pV1->getSeqLen());
            scaffGraph.addVertex(pScaffV1);
            numVerticedAdded++;
        }

        if (pScaffV2 == NULL)
        {
            pScaffV2 = new ScaffoldVertex(v2, pV2->getSeqLen());
            scaffGraph.addVertex(pScaffV2);
            numVerticedAdded++;
        }

        // Add an edge for this closure
        // Question: Is it okay for SGA Scaffolding algorithms if we have multiple edges between the same two scaffold vertexes?
        EdgeComp comp = (b.dir1 == b.dir2) ? EC_REVERSE : EC_SAME;
        int gapSize = iter->walks_[0].getEndToStartDistance();
        int gapStd = 0;
        ScaffoldLink link12(v2, b.dir1, comp, gapSize, gapStd, b.n, pV2->getSeqLen(), SLT_INFERRED);
        ScaffoldLink link21(v1, b.dir2, comp, gapSize, gapStd, b.n, pV1->getSeqLen(), SLT_INFERRED);
        ScaffoldEdge* pEdge12 = new ScaffoldEdge(pScaffV2, link12);
        ScaffoldEdge* pEdge21 = new ScaffoldEdge(pScaffV1, link21);
        pEdge12->setTwin(pEdge21);
        pEdge21->setTwin(pEdge12);
        scaffGraph.addEdge(pScaffV1, pEdge12);
        scaffGraph.addEdge(pScaffV2, pEdge21);
        numEdgesAdded += 2;
    }

    // Add edges to the scaffold graph
    cout << "Added " << numVerticedAdded << " vertices to scaffold graph.\n";
    cout << "Added " << numEdgesAdded << " edges to scaffold graph.\n";
   
    ScaffoldStatsVisitor statsVisit;
    scaffGraph.visit(statsVisit);

    std::cout << "Performing strict resolutions" << std::endl;

    // Remove transitive edges
    ScaffoldTransitiveReductionVisitor trVisit;
    scaffGraph.visit(trVisit);

    // Check for cycles in the graph
    ScaffoldAlgorithms::destroyStrictCycles(&scaffGraph, "scaffold.cycles.out");
    ScaffoldMultiEdgeRemoveVisitor meVisit;
    scaffGraph.visit(meVisit);
    scaffGraph.visit(statsVisit);

    // Linearize the scaffolds
    ScaffoldAlgorithms::makeScaffolds(&scaffGraph);

    // Break any remaining multi-edge contigs scaffolds
    ScaffoldMultiEdgeRemoveVisitor cutVisitor;
    scaffGraph.visit(cutVisitor);

    // Write out the scaffolds
    ScaffoldWriterVisitor writer("myScaffolds.txt");
    scaffGraph.visit(writer);

    ScaffoldCollectorVisitor scaffCollector;
    scaffGraph.visit(scaffCollector);
    vector<ScaffoldRecord> scaffs = scaffCollector.getScaffs();
    cout << "Have " << scaffs.size() << " scaffolds\n";

    // Collect each non-singleton scaffold and apply it to the graph
    vector<Closure> scaffClosures;
    typedef vector<ScaffoldRecord> ScaffRecVec;
    const ScaffRecVec::const_iterator iterB = scaffs.begin();
    const ScaffRecVec::const_iterator iterE = scaffs.end();
    for(ScaffRecVec::const_iterator iter = iterB; iter != iterE; iter++)
    {
        SGWalk walk = scaffoldToWalk(*iter);
        ostringstream scaffId;
        scaffId << "scaffold-" << (iter-iterB);
        scaffClosures.push_back(Closure(scaffId.str(), walk, 0, 0));
    }

    // Apply the closures to the graph
    cout << "Applying " << scaffClosures.size() << " scaffolds to the graph" << endl;
    ClosureAlgorithms::addClosuresToGraph(pGraph_, scaffClosures);
}

SGWalk ContigMerger::scaffoldToWalk(const ScaffoldRecord& rec)
{
    // Each edge in the scaffold graph corresponds to a walk in the 
    // string graph. Find the appropriate walk for each scaffold edge,
    // and orient the walk to match the orientation of the scaffold.
    string rootId = rec.getRoot();
    ScaffoldLinkVec links = rec.getLinks();
    EdgePtrVec scaffEdges;

//    cout << "-----------------" << endl;
//    cout << "KEYS:\n";
//    for (ClosureMap::const_iterator iter = closureMap_.begin();
//    iter != closureMap_.end();
//    iter++)
//    {
//        cout << iter->first << "\n";
//    }
//    cout << "------------------" << endl;

    for(size_t i = 0; i < links.size(); i++)
    {
        const ScaffoldLink& l = links[i];
        string nextId = l.endpointID;

        // Get the ClosureData for this ScaffoldLink
        ClosureDataKey key(rootId, l);
        //cout << "Looking for key: " << key << endl;
        ClosureMap::iterator iter = closureMap_.find(key);
        assert(iter != closureMap_.end());
        const ClosureData& data = iter->second;
        assert(data.walks_.size() == 1);

        const SGWalk& walk = data.walks_[0];

        // Get edges for the walk. Orient the edges so they are from 
        // the rootId to the nextId
        EdgePtrVec walkEdges = walk.getEdges();
        string id1 = walk.getStartVertex()->getID();
        string id2 = walk.getLastVertex()->getID();
        assert(walkEdges[0]->getStart()->getID() == id1);
        assert(walkEdges[walkEdges.size()-1]->getEnd()->getID() == id2);
        bool reverseWalk = (id1 == nextId);
        if (reverseWalk)
        {
            assert(id1 == nextId);
            assert(id2 == rootId);
        }
        else
        {
            assert(id1 == rootId);
            assert(id2 == nextId);
        }
        if (reverseWalk) walkEdges = reverse(walkEdges);

        // Copy walkEdges
        scaffEdges.insert(scaffEdges.end(), walkEdges.begin(), walkEdges.end());

        rootId = nextId;
    }

    // Check that the scaffold edges are sane
    if (scaffEdges.empty())
    {
       return SGWalk(pGraph_->getVertex(rec.getRoot()));
    }
    assert(scaffEdges.size() > 0);
    EdgeDir enterDir = scaffEdges[0]->getTwin()->getDir();
    const Vertex* nextV = scaffEdges[0]->getEnd();
    for (EdgePtrVec::const_iterator iter = scaffEdges.begin()+1;
         iter != scaffEdges.end();
         iter++)
    {
        const Edge* e = *iter;
        assert(e->getDir() == !enterDir);
        assert(nextV == e->getStart());
        enterDir = e->getTwin()->getDir();
        nextV = e->getEnd();
    }

    return SGWalk(scaffEdges);
}
