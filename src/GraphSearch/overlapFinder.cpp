#include <string>

#include "overlapFinder.h"
#include "bundle.h"
#include "OverlapTools.h"

using namespace std;

bool OverlapFinder::findOverlap(const Bundle* b, const StringGraph* pGraph, Overlap& overlap)
{
    // Calculate the range of allowed gaps
    float minGap = b->gap - maxStd_*b->std;
    float maxGap = b->gap + maxStd_*b->std;

    const float MIN_OL = 20.0;

    if (-minGap < MIN_OL)
        return false;

    // Convert the gaps to overlaps
    int minOverlap = -maxGap;
    if (minOverlap < MIN_OL) minOverlap = MIN_OL;
    int maxOverlap = -minGap;

    // Use the bundle to obtain the sequence. Note that the boundedOverlap function
    // computes overlaps of the suffix of s1 with the prefix of s2
    // s1: ----------------->
    // s2;            ------------->

    // There are four cases. To maintain sanity,
    // we will always use the sequence of vertex1 as s1.
    // (even if it means taking unnecessary reverse complements).
    string s1, s2, id1, id2;
    const Vertex* pVertex1 = pGraph->getVertex(b->vertex1ID);
    const Vertex* pVertex2 = pGraph->getVertex(b->vertex2ID);
    bool s1IsRC = false;
    bool s2IsRC = false;
    if ((b->dir1 == ED_SENSE) && (b->dir2 == ED_ANTISENSE))
    {
        s1 = pVertex1->getSeq().toString();
        s2 = pVertex2->getSeq().toString();
    }
    else if((b->dir1 == ED_SENSE) && (b->dir2 == ED_SENSE))
    {
        s1 = pVertex1->getSeq().toString();
        s2 = reverseComplementIUPAC(pVertex2->getSeq().toString());
        s2IsRC = true;
    }
    else if((b->dir1 == ED_ANTISENSE) && (b->dir2 == ED_SENSE))
    {
        s1 = reverseComplementIUPAC(pVertex2->getSeq().toString());
        s2 = reverseComplementIUPAC(pVertex1->getSeq().toString());
        s1IsRC = true;
        s2IsRC = true;
    }
    else if((b->dir1 == ED_ANTISENSE) && (b->dir2 == ED_ANTISENSE))
    {
        s1 = reverseComplementIUPAC(pVertex2->getSeq().toString());
        s2 = pVertex1->getSeq().toString();
        s1IsRC = true;
    }

    Match match;
    bool overlapFound = OverlapTools::boundedOverlapDP(s1, s2, minOverlap, maxOverlap, maxErrorRate_, match);

    if (!overlapFound)
        return false;

    // Restore the orientation of the match.  
    if (overlapFound)
    {
        assert(match.isReverse == false); 
        if (s1IsRC) match.coord[0].flip();
        if (s2IsRC) match.coord[1].flip();
        if (s1IsRC != s2IsRC) match.isReverse = true;
    }

    // Create the Overlap from the match:
    overlap = Overlap(b->vertex1ID, b->vertex2ID, match);
    return true;
}

