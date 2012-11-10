//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldRecord - A scaffold consisting of a 
// starting component and a vector of ordered
// links
//
#include "ScaffoldRecord.h"
#include "OverlapTools.h"
#include "SGSearch.h"

#define DEBUGRESOLVE 1
#define DEBUG_WALKTOPLACEMENT 1


void walkToContigPlacements(const SGWalk& walk, std::vector<ContigPlacement>& contigPlacements);

//
ScaffoldRecord::ScaffoldRecord() 
{

}

//
void ScaffoldRecord::setRoot(const std::string& root)
{
    m_rootID = root;
}

//
void ScaffoldRecord::addLink(const ScaffoldLink& link)
{
    m_links.push_back(link);
}

//
size_t ScaffoldRecord::getNumComponents() const
{
    if(m_rootID.empty())
        return 0;
    
    return 1 + m_links.size();
}

// Construct a string from the scaffold
std::string ScaffoldRecord::generateString(const ResolveParams& params, StringVector& contigPlacementDesc) const
{
    assert(params.pSequenceCollection != NULL);

    params.pStats->numScaffolds += 1;

    std::vector<ContigPlacement> contigPlacements;

    // Starting from the root, join the sequence(s) of the scaffold
    // together along with the appropriate gaps/overlap
    params.pSequenceCollection->setPlaced(m_rootID);
    std::string sequence = params.pSequenceCollection->getSequence(m_rootID);
    bool isRC = false;
    ContigPlacement rootPlacement(m_rootID, isRC, 0, sequence.size(), 0, sequence.size());
    contigPlacements.push_back(rootPlacement);

    if(m_links.empty())
    {
        contigPlacementDesc.push_back(rootPlacement.toString());
        return sequence;
    }

    EdgeDir rootDir = m_links[0].getDir();
    EdgeComp relativeComp = EC_SAME;
    EdgeComp prevComp = EC_SAME;
    std::string currID = m_rootID;

    // If this scaffold grows in the antisense direction,
    // we reverse every component and perform appends of the reversed
    // parts. After the scaffold is constructed we reverse again
    // to obtain the final scaffold in the desired orientation
    bool reverseAll = (rootDir == ED_ANTISENSE);
    if(reverseAll)
        sequence = reverse(sequence);

    // Iterate over all the linked contigs and append their sequence
    // to the scaffold
    for(size_t i = 0; i < m_links.size(); ++i)
    {   

        size_t currScaffPos = sequence.size();
        params.pStats->numGapsAttempted += 1;

        // Mark the current link as placed in the scaffold
        const ScaffoldLink& link = m_links[i];
        params.pSequenceCollection->setPlaced(link.endpointID);

        // Calculate the strand this sequence is on relative to the root
        if(link.getComp() == EC_REVERSE)
            relativeComp = !relativeComp;

        // Attempt to resolve the sequence of the link
        std::string resolvedSequence;
        

        // Step 1, try to walk through the graph between the vertices
        bool resolved = false;
        int insertedGap(0);
        if(params.resolveMask & RESOLVE_GRAPH_BEST || params.resolveMask & RESOLVE_GRAPH_UNIQUE)
        {
            SGWalk selectedWalk(NULL); // Create an empty walk with a NULL starting vertex.
            resolved = graphResolve(params, currID, link, resolvedSequence, selectedWalk);
            if(resolved)
            {

                // Determine the contig placements within the walk
                std::vector<ContigPlacement> newPlacements;
                walkToContigPlacements(selectedWalk, newPlacements);
                int walkL = selectedWalk.getStartToEndDistance();
                assert(newPlacements.back().scaffEnd_ == walkL);
   
                // Check that the contig placements either starts or ends with the current contig considered.
                // Walk will end with the current contig if the walk leaves the 5' end of the contig.
                bool contigIsFirst = (newPlacements[0].contigId_ == currID);
                bool contigIsLast = (newPlacements.back().contigId_ == currID);
                assert(contigIsFirst || contigIsLast);

                // Reverse the contig placements if necessary.
                if(contigIsLast)
                {
                    std::reverse(newPlacements.begin(), newPlacements.end());
                    for( size_t i = 0; i < newPlacements.size(); i++)
                        newPlacements[i].reverseScaffCoords(walkL);
                }

                // The returned sequence is wrt currID. If we flipped the sequence of currID, we must
                // flip this sequence
                if(prevComp == EC_REVERSE)
                {
                    resolvedSequence = reverseComplementIUPAC(resolvedSequence);
                    for( size_t i = 0; i < newPlacements.size(); i++)
                    {
                        newPlacements[i].reverseOrientation();
                    }
                }

                // Make sure the first placement is that of the current contig.
                // Then, remove this placement, and adjust all of the scaffold start and end coordinates
                assert(newPlacements[0].contigId_ == currID);
                int firstContigL = newPlacements[0].scaffEnd_;
                int adjustment = currScaffPos - firstContigL;
                for( size_t i = 1; i < newPlacements.size(); i++)
                {
                    ContigPlacement& placement = newPlacements[i];
                    placement.scaffStart_ += adjustment;
                    placement.scaffEnd_ += adjustment;
                }

                if(reverseAll)
                    resolvedSequence = reverse(resolvedSequence);

                params.pStats->numGapsResolved += 1;

                // Save new contig placements. Note that the first placement is that of a contig
                // we hace placed previously
                assert(newPlacements[0].contigId_ == currID);
                assert(contigPlacements.back().contigId_ == currID);
                contigPlacements.insert(contigPlacements.end(), newPlacements.begin()+1, newPlacements.end());
            }
        }

        // Step 2, try to resolve a predicted overlap between current sequence
        // and the sequence of the linked contig
        if(!resolved)
        {
            currID = link.endpointID;
            bool isRC = (relativeComp == EC_REVERSE);
            // Get the sequence that should be potentially appended in
            std::string contigSeq = params.pSequenceCollection->getSequence(currID);
            size_t contigSeqL = contigSeq.size();

            if(relativeComp == EC_REVERSE)
                contigSeq = reverseComplementIUPAC(contigSeq);
            if(reverseAll)
                contigSeq = reverse(contigSeq);
            
            //
            if(link.distance < 0 && params.resolveMask & RESOLVE_OVERLAP)
            {
                resolved = overlapResolve(params, sequence, contigSeq, link, resolvedSequence);
                if(resolved)
                {
                    params.pStats->numGapsResolved += 1;
                    params.pStats->overlapFound += 1;

                    // Add the placement of the contig
                    size_t resolvedSeqL = resolvedSequence.size();
                    size_t scaffEnd = currScaffPos + resolvedSeqL;
                    int contigSeqEnd =  contigSeqL;
                    int contigSeqStart = contigSeqEnd - resolvedSeqL;
                    // If contig is reverse-complement, select the contig seq start/end coords
                    // with respect to the contig's forward orientation
                    if (isRC)
                    {
                        contigSeqStart = 0;
                        contigSeqEnd = resolvedSeqL;
                    }
                    contigPlacements.push_back(ContigPlacement(currID, isRC, currScaffPos, scaffEnd, contigSeqStart, contigSeqEnd));
                }
                else
                {
                    params.pStats->overlapFailed += 1;
                }
            }

            // Step 3, just introduce a gap between the sequences
            if(!resolved)
            {
                introduceGap(params.minGapLength, contigSeq, link, resolvedSequence, insertedGap);
                // If resolved by gap, adjust the seqStart to be the end of the gap
                assert(insertedGap > 0);
                size_t scaffStart = currScaffPos + insertedGap; // Start position of contig in scaffold, after the gap
                size_t scaffEnd = scaffStart + resolvedSequence.size();
                contigPlacements.push_back(ContigPlacement(currID, isRC, scaffStart, scaffEnd, 0, contigSeqL));
            }
        }
        sequence.append(resolvedSequence);
        prevComp = relativeComp;
    }

    if(reverseAll) 
    {
        sequence = reverse(sequence);
        size_t scaffLength = sequence.size();
        std::reverse(contigPlacements.begin(), contigPlacements.end());
        for (size_t i = 0; i < contigPlacements.size(); i++)
        {
            contigPlacements[i].reverseScaffCoords(scaffLength);
        }
    }

    contigPlacementDesc.clear();
    for (size_t i = 0; i < contigPlacements.size(); i++)
        contigPlacementDesc.push_back(contigPlacements[i].toString());

    return sequence;
}

// Attempt to resolve a scaffold link by finding a walk through the graph linking the two vertices
bool ScaffoldRecord::graphResolve(const ResolveParams& params, const std::string& startID, 
                                  const ScaffoldLink& link, std::string& outExtensionString,
                                  SGWalk& selectedWalk) const
{
    assert(params.pGraph != NULL);

    // Get the vertex to start the search from
    Vertex* pStartVertex = params.pGraph->getVertex(startID);
    Vertex* pEndVertex = params.pGraph->getVertex(link.endpointID);
    assert(pStartVertex != NULL && pEndVertex != NULL);

    int threshold = static_cast<int>(params.distanceFactor * link.stdDev);
    int maxDistance = link.distance + threshold;
    int maxExtensionDistance = maxDistance + pEndVertex->getSeqLen();
    SGWalkVector walks;
    SGSearch::findWalks(pStartVertex, pEndVertex, link.getDir(), maxExtensionDistance, 10000, true, walks);

    int numWalksValid = 0;
    int numWalksClosest = 0;
    int selectedIdx = -1;
    int closestDist = std::numeric_limits<int>::max();

#ifdef DEBUGRESOLVE
            std::cout << "Attempting graph resolve of link " << startID << " -- " << link.endpointID << " expected distance: " << link.distance << " orientation: " << link.edgeData.getComp() << "\n";
#endif
    
    // Select the closest walk to the distance estimate
    for(size_t i = 0; i < walks.size(); ++i)
    {
        // Check that the orientation of the walk is the same as the expected link
        std::vector<EdgeComp> vertexOrientations = walks[i].getOrientationsToStart();
        assert(walks[i].getLastEdge()->getEndID() == link.endpointID);

        if(vertexOrientations.back() != link.edgeData.getComp())
        {
#ifdef DEBUGRESOLVE
            std::cout << "SKIPPING WALK OF THE WRONG ORIENTATION\n";
#endif
            continue;
        }

        int walkDistance = walks[i].getEndToStartDistance();
        int diff = abs(abs(link.distance - walkDistance));
        if(diff <= threshold)
        {

#ifdef DEBUGRESOLVE
            std::cout << "  Walk distance: " << walkDistance << " diff: " << diff << " threshold: " << threshold << " close: " << closestDist << "\n";
#endif
            ++numWalksValid;
            if(diff < closestDist)
            {
                selectedIdx = i;
                closestDist = diff;
                numWalksClosest = 1;
            }
            else if(diff == closestDist)
            {
                numWalksClosest += 1;
            }

        }
    }

    // Choose the best path, if any, depending on the algorithm to use
    bool useWalk = false;

    if(numWalksValid > 0)
    {
        if(params.resolveMask & RESOLVE_GRAPH_BEST)
        {
            // If the unique flag is not set, or we only have 1 closest walk, select it
            if(!(params.resolveMask & RESOLVE_GRAPH_UNIQUE) || numWalksClosest == 1)
                useWalk = true;
            else if((params.resolveMask & RESOLVE_GRAPH_UNIQUE) && numWalksClosest > 1)
                params.pStats->graphWalkTooMany += 1;
        }
        else
        {
            if(numWalksValid == 1)
                useWalk = true;
            else if(numWalksValid > 1)
                params.pStats->graphWalkTooMany += 1;
        }
    }

#ifdef DEBUGRESOLVE    
    std::cout << "  Num walks: " << walks.size() << " Num valid: " << numWalksValid << " Num closest: " << numWalksClosest << " using: " << useWalk << "\n";
#endif

    // Was an acceptable walk found? 
    if(useWalk)
    {
        assert(selectedIdx != -1);
        selectedWalk = walks[selectedIdx];
        outExtensionString = selectedWalk.getString(SGWT_EXTENSION);
        params.pStats->graphWalkFound += 1;

        // Mark all vertices in the walk as visited
        VertexPtrVec vertexPtrVector = selectedWalk.getVertices();
        for(size_t i = 0; i < vertexPtrVector.size(); ++i)
            params.pSequenceCollection->setPlaced(vertexPtrVector[i]->getID());
        return true;
    }
    else
    {
        if(numWalksValid == 0)
            params.pStats->graphWalkNoPath += 1;
        assert(outExtensionString.empty());
        return false;
    }
}

// Attempt to resolve a predicted overlap between s1 and s2
// Returns true if there overlap was found and the overhang of s2 is placed in outString
bool ScaffoldRecord::overlapResolve(const ResolveParams& params, const std::string& s1, const std::string& s2, 
                                    const ScaffoldLink& link, std::string& outString) const
{
    // Attempt to find an overlap between these sequences
    int expectedOverlap = -1 * link.distance;

#ifdef DEBUGRESOLVE
    std::cout << "Attempting overlap resolve of link to " << link.endpointID << " expected distance: " << link.distance << " orientation: " << link.edgeData.getComp() << "\n";
#endif


    // If the maximum overlap was not set, set it to the expected overlap * 3 stddev
    int upperBound = 0;
    if(params.maxOverlap == -1)
        upperBound = static_cast<int>(expectedOverlap + 3.0f * link.stdDev);
    else
        upperBound = params.maxOverlap;
    
    // Calculate the best match
    Match match;
    bool overlapFound = OverlapTools::boundedOverlapDP(s1, s2, params.minOverlap, upperBound, params.maxErrorRate, match);
    if(overlapFound)
    {
#ifdef DEBUGRESOLVE
        std::cout << "Overlap found, length: " << match.coord[1].length() << "\n";
#endif
        SeqCoord overlapCoord = match.coord[1];
        SeqCoord overhangCoord = overlapCoord.complement();
        outString = overhangCoord.getSubstring(s2);
        return true;
    }
    else
    {
        return false;
    }
}

// Resolve a link with a gap
bool ScaffoldRecord::introduceGap(int minGapLength, const std::string& contigString, const ScaffoldLink& link, std::string& out,
                                  int& insertedGap) const
{
    assert(out.empty());
    if(link.distance < 0)
    {
        // Truncate the string using the expected overlap and add a gap with a fixed number of Ns
        out.append(minGapLength, 'N');
        insertedGap = minGapLength;
        int expectedOverlap = -1 * link.distance;
        assert(expectedOverlap < (int)contigString.length());
        out.append(contigString.substr(expectedOverlap));
    }
    else
    {
        int gap = std::max(link.distance, minGapLength);
        insertedGap = gap;
        out.append(gap, 'N');
        out.append(contigString);
    }
    return true;
}

//
void ScaffoldRecord::parse(const std::string& text)
{
    StringVector fields = split(text, '\t');
    assert(fields.size() >= 1);

    m_rootID = fields[0];
    for(size_t i = 1; i < fields.size(); ++i)
    {
        ScaffoldLink link;
        link.parse(fields[i]);
        m_links.push_back(link);
    }
}

//
void ScaffoldRecord::writeScaf(std::ostream* pWriter)
{
    *pWriter << m_rootID;
    for(size_t i = 0; i < m_links.size(); ++i)
        *pWriter << "\t" << m_links[i];
    *pWriter << "\n";
}

/////////////////////////////////////////////////////////////////////////////////
// Given a walk resolving a scaffold gap, determine the contig placement information
// for all contigs of the walk.
// Position 0 is relative to the beginning of the walk

// Example: walk from A to B:
//   A |------------->                         <---------------| B
//                |-------------> C
//                         <----------| D
//                                 <-----------------| E
//                                          |------------> F
// Placements:
//                   |--C--------|-D---|----E---|---------------| B
// Note:
// We do not plact contig A
// Both D & E are reverse-complement.
// F does not get placed because we use all of B.
void walkToContigPlacements(const SGWalk& walk, std::vector<ContigPlacement>& contigPlacements)
{
    contigPlacements.clear();

    // Get the sequence of the walk and the placements of all contigs in the walk
    // Position 0 is the start of the first contig in the walk
    SGWalkVertexPlacementVector vertexPlacements;
    std::string walkSeq = walk.getString(SGWT_START_TO_END, &vertexPlacements);

    // Add the placement of the first contig of walk
    SGWalkVertexPlacement& vPlacement = vertexPlacements[0];
    Vertex * pFirstVertex = vPlacement.pVertex;
    size_t firstContigL = pFirstVertex->getSeqLen();
    assert(vPlacement.position == 0);
    ContigPlacement contigPlacement(pFirstVertex->getID(), vPlacement.isRC, 0, firstContigL, 0, firstContigL);
    contigPlacements.push_back(contigPlacement);

    int currPos = firstContigL; // absolute position in walk
    int gapStart = currPos;
    int gapSize = walk.getEndToStartDistance();
    int gapEnd = gapStart + (gapSize > 0 ? gapSize : 0); // absolute position in walk

    // Add placements from interior contigs if they contribute to the walk's sequence
    // by filling in a gap between the first and last contig
    size_t lastContig = vertexPlacements.size()-1;
    for (size_t i = 1; i < lastContig; i++)
    {
        if (currPos >= gapEnd) break;

        SGWalkVertexPlacement& vPlacement = vertexPlacements[i];
        Vertex * pVertex = vPlacement.pVertex;

        // Get the contig start/ending position in walk
        int contigStartPos = vPlacement.position;
        int contigSeqL = pVertex->getSeqLen();
        int contigEndPos = contigStartPos + contigSeqL;

        // Coordinates of walk interval first covered by this contig
        int scaffSeqStart = currPos;
        int scaffSeqEnd = (contigEndPos > gapEnd) ? gapEnd : contigEndPos;
        int seqContributionL = scaffSeqEnd - scaffSeqStart;

        // Coordinates of contig matched to this walk interval
        int contigSeqStart = currPos - contigStartPos;
        int contigSeqEnd = contigSeqStart + seqContributionL;
        assert(contigSeqStart > 0);
        assert(contigSeqEnd <= contigSeqL);

        // If contig is reverse-complement, select the contig seq start/end coords
        // with respect to the contig's forward orientation
        if (vPlacement.isRC)
        {
            contigSeqStart = contigSeqL - contigSeqEnd;
            contigSeqEnd = contigSeqStart + seqContributionL;
        }
        ContigPlacement contigPlacement(pVertex->getID(), vPlacement.isRC, scaffSeqStart, scaffSeqEnd,
                                        contigSeqStart, contigSeqEnd);
        contigPlacements.push_back(contigPlacement);
        currPos += seqContributionL;
    }

    //Add placement of final contig
    assert(currPos == gapEnd);
    vPlacement = vertexPlacements[lastContig];
    Vertex * pLastVertex = vPlacement.pVertex;
    // Get the contig start/ending position in walk
    int contigStartPos = vPlacement.position;
    int contigSeqL = pLastVertex->getSeqLen();
    int contigEndPos = contigStartPos + contigSeqL;
    // Coordinate of contig matched to scaffold sequence
    int contigSeqStart = currPos - contigStartPos;
    int contigSeqEnd = contigSeqL;
    // Coordinate of scaffold matched to contig sequence
    int scaffSeqStart = currPos;
    int scaffSeqEnd = contigEndPos;
    contigPlacement = ContigPlacement(pLastVertex->getID(), vPlacement.isRC, scaffSeqStart, scaffSeqEnd, contigSeqStart, contigSeqEnd);
    contigPlacements.push_back(contigPlacement);
    return;
}
