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

//#define DEBUGRESOLVE 1
//#define DEBUGPLACEMENTS 1


void walkToContigPlacements(const SGWalk& walk, const ContigPlacement& currPlacement, const std::string& nextID,
                            std::vector<ContigPlacement>& newPlacements);

void printPlacements(const std::vector<ContigPlacement>& placements)
{
    for (size_t i = 0; i < placements.size(); i++)
        std::cout << "\t" << placements[i].toString() << std::endl;
}

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

        const size_t currScaffPos = sequence.size();
        params.pStats->numGapsAttempted += 1;

        // Mark the current link as placed in the scaffold
        const ScaffoldLink& link = m_links[i];
        const std::string nextID = link.endpointID; // The ID of the next sequence in scaffold we are trying to place
        params.pSequenceCollection->setPlaced(nextID);

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
                assert(currID == contigPlacements.back().contigId_);
                const ContigPlacement& currPlacement = contigPlacements.back();
                walkToContigPlacements(selectedWalk, currPlacement, nextID, newPlacements);


                #if DEBUGPLACEMENTS
                std::cout << "GAP RESOLVED WITH GRAPH WALK."
                          << " CURRID=" << currID
                          << " NEXTID=" << nextID
                          << "\nNEW CONTIG PLACEMENTS:\n";
                printPlacements(newPlacements);
                #endif

                // The returned sequence is wrt currID. If we flipped the sequence of currID, we must
                // flip this sequence
                if(prevComp == EC_REVERSE)
                {
                    resolvedSequence = reverseComplementIUPAC(resolvedSequence);
                }

                if(reverseAll)
                    resolvedSequence = reverse(resolvedSequence);

                contigPlacements.insert(contigPlacements.end(), newPlacements.begin(), newPlacements.end());

                params.pStats->numGapsResolved += 1;
            }
        }

        // Step 2, try to resolve a predicted overlap between current sequence
        // and the sequence of the linked contig
        if(!resolved)
        {
            bool isRC = (relativeComp == EC_REVERSE);
            // Get the sequence that should be potentially appended in
            std::string contigSeq = params.pSequenceCollection->getSequence(nextID);
            const size_t contigSeqL = contigSeq.size();

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
                    size_t resolvedSeqL = resolvedSequence.size();
                    size_t scaffEnd = currScaffPos + resolvedSeqL;

                    // Check which end of contig ends up in the overlap.
                    // Case 1: !isRC && !reverseAll:
                    // scaff start -------> ....
                    //                            ------> contig (overlap at front)
                    // Case 2: !isRC && reverseAll:
                    // contig -----> (overlap at end)
                    //                .....  ------> scaff start 
                    // Case 3: isRC && !reverseAll:
                    // ------> scaff start
                    //         ..... <----- contig (overlap at end)
                    // Case 4: isRC && reverseAll:
                    //  <----- contig (overlap at front)
                    //        ....  -------> scaff start
                    bool overlapInFront = (!isRC && !reverseAll) || (isRC && reverseAll);
                    int contigSeqStart, contigSeqEnd;
                    if (overlapInFront)
                    {
                        contigSeqEnd =  contigSeqL;
                        contigSeqStart = contigSeqEnd - resolvedSeqL;
                    }
                    else
                    {
                        contigSeqStart = 0;
                        contigSeqEnd = resolvedSeqL;
                    }

                    ContigPlacement newPlacement(nextID, isRC, currScaffPos, scaffEnd, contigSeqStart, contigSeqEnd);
                    contigPlacements.push_back(newPlacement);
                    #if DEBUGPLACEMENTS
                    std::cout << "RESOLVED LINK WITH OVERLAP. NEW PLACEMENT:"
                              << "\t" << newPlacement.toString() << std::endl;
                    #endif
                }
                else
                {
                    params.pStats->overlapFailed += 1;
                }
            }

            // Step 3, just introduce a gap between the sequences
            if(!resolved)
            {

                // Reminder of how introduceGap works:
                // - If there is a predicted overlap, then minGapLength N's are added, and then
                //   contig sequence is trimmed by this number of bases before insertion
                // - If there is not a predicted overlap, then min(minGapLength, gapLength) N's are
                //   added, followed by the FULL sequence of contig
                introduceGap(params.minGapLength, contigSeq, link, resolvedSequence, insertedGap);

                // If resolved by gap, adjust the seqStart to be the end of the gap
                assert(insertedGap > 0);

                // Add the placement of contig nextID
                size_t scaffStart = currScaffPos + insertedGap; // Start position of contig in scaffold, after the gap
                size_t scaffEnd = currScaffPos + resolvedSequence.size();
                size_t numContigBasesUsed = scaffEnd - scaffStart;
                assert(numContigBasesUsed <= contigSeqL);

                // Check which end of contig ends up in the gap. If introduceGap() clips sequence
                // off of the contig, we must carefully record this
                // Case 1: !isRC && !reverseAll:
                // scaff start -------> ....
                //                            ------> contig (gap at front)
                // Case 2: !isRC && reverseAll:
                // contig -----> (gap at end)
                //                .....  ------> scaff start 
                // Case 3: isRC && !reverseAll:
                // ------> scaff start
                //         ..... <----- contig (gap at end)
                // Case 4: isRC && reverseAll:
                //  <----- contig (gap at front)
                //        ....  -------> scaff start
                bool clipInFront = (!isRC && !reverseAll) || (isRC && reverseAll);
                int contigSeqStart, contigSeqEnd;
                if (clipInFront)
                {
                    contigSeqEnd =  contigSeqL;
                    contigSeqStart = contigSeqEnd - numContigBasesUsed;
                }
                else
                {
                    contigSeqStart = 0;
                    contigSeqEnd = numContigBasesUsed;
                }

                ContigPlacement newPlacement(nextID, isRC, scaffStart, scaffEnd, contigSeqStart, contigSeqEnd);
                contigPlacements.push_back(newPlacement);
                #if DEBUGPLACEMENTS
                std::cout << "FAILED TO RESOLVE WITH OVERLAP. PLACEMENT AFTER GAP:"
                          << "\t" << newPlacement.toString() << std::endl;
                #endif
            }
        }

        sequence.append(resolvedSequence);
        prevComp = relativeComp;
        currID = nextID;
    }

    #if DEBUGPLACEMENTS
    std::cout << "DONE MAKING ALL CONTIG PLACEMENTS IN SCAFFOLD." << std::endl;
    #endif

    if(reverseAll) 
    {
        sequence = reverse(sequence);
        size_t scaffLength = sequence.size();

        #if DEBUGPLACEMENTS
        std::cout << "MUST REVERSE PLACEMENTS OF CONTIGS IN SCAFFOLDS.\n"
                  << "ScaffLength: " << scaffLength << "\n"
                  << "PLACEMENTS BEFORE REVERSING:" << std::endl;
        printPlacements(contigPlacements);
        #endif

        std::reverse(contigPlacements.begin(), contigPlacements.end());
        for (size_t i = 0; i < contigPlacements.size(); i++)
            contigPlacements[i].reverseScaffCoords(scaffLength);

        #if DEBUGRESOLVE
        std::cout << "PLACEMENTS AFTER REVERSING: \n";
        printPlacements(contigPlacements);
        #endif
    }


    // Form the contig placement description string vector
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
// This is a helper function for ScaffoldRecord::generateString().
//
// Given a walk resolving a scaffold gap from currID to nextID,
// determine the contig placement information for all contigs of the walk.
// The contigs are placed in a non-overlapping fashion to
// completely cover the path, giving priority first to currID, second to nextID,
// and then to all other contigs in the gap, in order of walk from currID to nextID.
//
// Orient placed contigs so that they are consistent relative to the orientation of
// contig currID in the scaffold (as specified by currPlacement)
//
// Adjust all contig placement positions so that they extend the placement of contig
// currID in the scaffold (as specified in currPlacement)

// Example: currID = A, nextID = B
//   A |------------->                         <---------------| B
//                |-------------> C
//                         <----------| D
//                                 <-----------------| E
//                                          |------------> F
// Placements:
//     |------A-------|----C----|-D---|----E---|---------------| B
//     0
//
// Note:
// We place contig C starting with the first base beyond the end of contig A, since A is currID and gets priority.
// F does not get placed because we give priority to B (which is nextID), and F does not contribute any novel sequence that isn't covered by B or preceeding contigs.
//
// Another Example: currID = A, nextID = B
//   B |------------->      
//                 |---------------> A
//             |-------------> C
//
// Placements:
//     |------B----|--------A------|
//                                 0
// Note:
// A & B overlap directly, so we do not place C.
// A is currID, so it gets priority over B.
void walkToContigPlacements(const SGWalk& walk, const ContigPlacement& currPlacement, const std::string& nextID,
                            std::vector<ContigPlacement>& newPlacements)
{

    std::vector<ContigPlacement> contigPlacements;

    std::string currID = currPlacement.contigId_;

    // Get the sequence of the walk and the placements of all contigs in the walk
    // Position 0 is the start of the first contig in the walk
    SGWalkVertexPlacementVector vertexPlacements;
    std::string walkSeq = walk.getString(SGWT_START_TO_END, &vertexPlacements);
    int walkL = walk.getStartToEndDistance();
    assert(walkL == (int) walkSeq.size());
    #if DEBUGPLACEMENTS
    std::cout << "walkToContigPlacements:\n"
              << "\tcurrID=" << currID << "\n"
              << "\tnextID=" << nextID << "\n"
              << "\twalkLength=" << walkL << "\n"
              << "\twalk=:\n";
    walk.print();
    #endif

    int currContigL;
    bool currIsFirst = (vertexPlacements.front().pVertex->getID() == currID);
    bool currIsLast = (vertexPlacements.back().pVertex->getID() == currID);
    assert(currIsFirst || currIsLast);
    if (currIsFirst)
    {
        currContigL = vertexPlacements.front().pVertex->getSeqLen();
        assert(nextID == vertexPlacements.back().pVertex->getID());
    }
    if (currIsLast)
    {
        currContigL = vertexPlacements.back().pVertex->getSeqLen();
        assert(nextID == vertexPlacements.front().pVertex->getID());
    }
    
    // If currID is the last contig placed in the walk, take the reverse complement of the
    // walk such that we start with current contig at position 0
    if (currIsLast)
    {
        std::reverse(vertexPlacements.begin(), vertexPlacements.end());
        for (size_t i = 0; i < vertexPlacements.size(); i++)
        {
            SGWalkVertexPlacement& vp = vertexPlacements[i];
            vp.isRC = !vp.isRC;
            int newPos = walkL - (vp.position + vp.pVertex->getSeqLen());
            assert(newPos >= 0);
            assert(newPos < walkL);
            vp.position = newPos;
        }
    }
    assert(vertexPlacements.front().pVertex->getID() == currID);
    assert(vertexPlacements.back().pVertex->getID() == nextID);

    /////////////////////////////////////////////////////////////////
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

    /////////////////////////////////////////////////////////////////
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
        assert(scaffSeqStart > 0);
        assert(scaffSeqEnd > scaffSeqStart);
        assert(seqContributionL > 0);

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

    ////////////////////////////////////////////////////////////
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
    int seqContributionL = scaffSeqEnd - scaffSeqStart;
    // If contig is reverse-complement, select the contig seq start/end coords
    // with respect to the contig's forward orientation
    if (vPlacement.isRC)
    {
        contigSeqStart = contigSeqL - contigSeqEnd;
        contigSeqEnd = contigSeqStart + seqContributionL;
    }
    contigPlacement = ContigPlacement(pLastVertex->getID(), vPlacement.isRC, scaffSeqStart, scaffSeqEnd, contigSeqStart, contigSeqEnd);
    contigPlacements.push_back(contigPlacement);


    #if DEBUGPLACEMENTS
    std::cout << "CONTIG PLACEMENTS IN WALK FROM CURRID="
              << currID << " TO NEXTID=" << nextID
              << " (BEFORE ADJUSTING POSITIONS):\n";
    printPlacements(contigPlacements);
    #endif

    ////////////////////////////////////////////////////////
    // ASSERT SANITY OF THE PLACEMENTS
    assert(currID == contigPlacements.front().contigId_);
    assert(nextID == contigPlacements.back().contigId_);
    assert(0 == contigPlacements.front().scaffStart_);
    assert(contigPlacements.back().scaffEnd_ == walkL);

    // Check that the placements are ungapped
    int lastEnd = contigPlacements.front().scaffEnd_;
    for (size_t i = 1; i < contigPlacements.size(); i++)
    {
        const ContigPlacement& p = contigPlacements[i];
        assert(p.scaffStart_ == lastEnd);
        int contigSpan = p.contigEnd_ - p.contigStart_;
        int scaffSpan = p.scaffEnd_ - p.scaffStart_;
        assert(contigSpan == scaffSpan);
        assert(contigSpan > 0);
        assert(scaffSpan > 0);
        lastEnd = p.scaffEnd_;
    }
    /////////////////////////////////////////////////////////////

    // Adjust the scaffold start and end positions relative to the placement
    // of contig currID in the full scaffold, as given by currPlacement
    // If the orientation of contig currID in contigPlacements does not match
    // the orientation in the full scaffold (as given by currPlacement),
    // reverse all orientations in contigPlacements
    assert(currID == contigPlacements[0].contigId_);
    bool reverseOrientations = (contigPlacements[0].isRC_ != currPlacement.isRC_);
    int currScaffPos = currPlacement.scaffEnd_;
    int adjustment = currScaffPos - currContigL;
    for( size_t i = 0; i < contigPlacements.size(); i++)
    {
        ContigPlacement& placement = contigPlacements[i];
        placement.scaffStart_ += adjustment;
        placement.scaffEnd_ += adjustment;
        if (reverseOrientations)
            placement.reverseOrientation();
    }

    // Copy all of the contig placements (excluding placement for currID)
    // into new placements
    newPlacements.clear();
    newPlacements.insert(newPlacements.end(), contigPlacements.begin()+1, contigPlacements.end());
}
