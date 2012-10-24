#  File: investigateSam.py
#  Author: Lee Mendelowitz
#  Date: 8/1/2012
#
#  Tools for iterating over a SAM/BAM file and computing
#  contig coverage, insert size histogram, and generating links
#  between contigs for scaffolding
###############################################################

import sys
import pysam
import numpy as np
from histogram import Histogram
from collections import defaultdict

###############################################
# Yield pairs from an open pysam.SamFile
# This assumes that the alignment records in the samFile sorted in pairs.
# Ex:
# [alignment for insert0, read 1]
# [alignment for insert0, read 2]
# [alignment for insert1, read 1]
# [alignment for insert1, read 2]
# etc.
def yieldPairs(samfile):
    try:
        while True:
            read1,read2 = samfile.next(), samfile.next()
            (r1,r2) = (read1,read2) if read1.is_read1 else (read2,read1)
            yield (r1,r2)
    except StopIteration:
        pass


###############################################
# Get the base name for a read in an alignment record
read_exts = set(['.1', '.2', '/1', '/2'])
def getReadBn(readName):
    if readName[-2:] in read_exts:
        return(readName[:-2])
    return readName

###############################################
# Yield pairs from an open pysam.SamFile
# The alignments should be sorted in order of read name
# The assumes that the read names:
#   a) end in '/1' and '/2'
#   b) or end in '.1' and '.2'
#   c) or are identical
# This will skip any reads that are aligned as singletons (were read pair is missing alignment)
def yieldPairsSmart(samFile):
    read1 = None
    read2 = None
    try:
        while True:
            if not read1: read1 = samFile.next()
            if not read2: read2 = samFile.next()
            if (getReadBn(read1.qname) == getReadBn(read2.qname)) and \
               (not read1.is_unmapped) and \
               (not read2.is_unmapped):
               yield (read1, read2)
               read1, read2 = None, None
            else:
                read1, read2 = (None, None) if read2.is_unmapped else (read2, None)
    except StopIteration:
        pass

###############################################
def pairIsValidSameContig(r1,r2):
    if r1.is_unmapped or r2.is_unmapped:
        return False
    
###############################################
def getInsSize(samFileName, fileType='sam'):
    insSize = []
    invalidCount = 0
    pairCount = 0
    diffContigs = 0

    if fileType == 'sam':
        samFile = pysam.Samfile(samFileName, 'rb')
    else:
        samFile = pysam.Samfile(samFileName, 'r')

    for r1, r2 in yieldPairs(samFile):
        pairCount += 1

        if r1.is_unmapped or r2.is_unmapped:
            invalidCount += 1
            continue

        if r1.tid != r2.tid:
            diffContigs += 1
            continue

        # Being here means that the reads are mapped to the same contig
        if r1.pos < r2.pos:
            (leftRead, rightRead) = (r1, r2)
        elif r2.pos < r1.pos:
            (leftRead, rightRead) = (r2, r1)
        else:
            (leftRead, rightRead) = (r1, r2) if not (r1.is_reverse) else (r2, r1)

        if leftRead.is_reverse or not rightRead.is_reverse:
            invalidCount += 1
            continue
        else:
            insSize.append(abs(leftRead.tlen))
    print 'Read %i pairs.\n%i pairs are invalid.\n%i pairs are in different contigs.'%(pairCount, invalidCount, diffContigs)
    return insSize


###############################################
# Read BAM File. Generate histogram of insert size and links between contigs
def readBam(bamFileName, insLow, insHigh, insInc, fileType='bam'):

    if fileType == 'bam':
        bamFile = pysam.Samfile(bamFileName, 'rb')
    else:
        bamFile = pysam.Samfile(bamFileName, 'r')

    contigLengthDict = dict(zip(bamFile.references, bamFile.lengths))

    binEdges = np.arange(insLow, insHigh+1, insInc)
    H =  Histogram(binEdges)
    invalidCount = 0
    pairCount = 0
    diffContigs = 0
    links = [] # Links between different contigs
    for r1, r2 in yieldPairs(bamFile):
        pairCount += 1

        if r1.is_unmapped or r2.is_unmapped:
            invalidCount += 1
            continue

        if r1.tid != r2.tid:
            diffContigs += 1
            # Create a link between contigs
            # Assume read1 is forward and read2 is reverse
            # Therefore, if read1 is mapped as forward to contig1, then contig1 is forward
            # and if read2 is mapped as forward to contig2, then contig2 is reverse.
            contig1 = bamFile.getrname(r1.tid)
            contig1L = contigLengthDict[contig1]
            contig2 = bamFile.getrname(r2.tid)
            contig2L = contigLengthDict[contig2]
            or1 = 'F' if not r1.is_reverse else 'R'
            or2 = 'F' if r2.is_reverse else 'R'

            # d1 and d2 capture the amount of the insert size that are contained in each contig.
            # |-------------------|
            #       r1 |--->
            #          |<----d1-->|
            if not r1.is_reverse:
                or1 = 'F'
                d1 = int(contig1L - r1.pos) 
            else:
                or1 = 'R'
                d1 = int(r1.positions[-1])+1

            if not r2.is_reverse:
                or2 = 'R'
                d2 = int(contig2L - r2.pos)
            else:
                or2 = 'F'
                d2 = int(r2.positions[-1])+1
            links.append((contig1, or1, d1, contig2, or2, d2))
            continue

        # Being here means that the reads are mapped to the same contig
        if r1.pos < r2.pos:
            (leftRead, rightRead) = (r1, r2)
        elif r2.pos < r1.pos:
            (leftRead, rightRead) = (r2, r1)
        else:
            (leftRead, rightRead) = (r1, r2) if not (r1.is_reverse) else (r2, r1)

        if leftRead.is_reverse or not rightRead.is_reverse:
            invalidCount += 1
            continue
        else:
            H.add(abs(leftRead.tlen))
    print 'Read %i pairs.\n%i pairs are invalid.\n%i pairs are in different contigs.'%(pairCount, invalidCount, diffContigs)
    return H, links

###############################################
# Read BAM File. Generate histogram of insert size.
# Return a dictionary of links between contigs. 
# The linkDict return has keys: (contig1, contig2) where contig1 <= contig2
#                       values: list of (or1, d1, or2, d2) tuples where d is the amount of
#                               insert contained in each contig
#
#                          Example of d for read1 if read1 orientation is supposed to be forward:
#                          |---------------> contig
#                                   |--> read1
#                                   |<---->| d
# 
def readBam_ver2(bamFileName, insLow, insHigh, insInc, r1Or = 'F', r2Or = 'R', fileType='bam'):
    assert(r1Or in 'FR' and len(r1Or)==1)
    assert(r2Or in 'FR' and len(r2Or)==1)
    r1IsForward = (r1Or == 'F')
    r2IsForward = (r2Or == 'F')

    if fileType == 'bam':
        bamFile = pysam.Samfile(bamFileName, 'rb')
    else:
        bamFile = pysam.Samfile(bamFileName, 'r')

    contigLengthDict = dict(zip(bamFile.references, bamFile.lengths))

    # Set up Histogram
    binEdges = np.arange(insLow, insHigh+1, insInc)
    H =  Histogram(binEdges)

    pairCount = 0
    linkCount = 0
    concordantCount = 0

    links = [] # Links between different contigs

    linkDict = defaultdict(list)

    #############################################################################
    # Set up functions to compute the orientation of contig and the
    # length of insert contained:
    # Use readSense() to measure insert and orientation based on how read points:
    # |--------------------> contig
    #               |-->  read
    #               |<---->| d
    # OR
    # |-------------------> contig
    #     |<--  read
    # |<---->| d
    def readSense(rec):
        if rec.is_reverse:
            # Read is reverse wrt contig
            return ('B', int(rec.positions[-1])+1)
        else: 
            # Read is forward wrt contig
            cL = contigLengthDict[bamFile.getrname(rec.tid)]
            return ('E', int(cL - rec.pos))

    # Use readAntisense() to measure insert and orientation based on the opposite of how read points:
    # |--------------------> contig
    #               <--|  read
    #               |<---->| d
    # OR
    # |-------------------> contig
    #     |-->  read
    # |<---->| d
    def readAntisense(rec):
        if rec.is_reverse:
            # Read is reverse wrt contig
            cL = contigLengthDict[bamFile.getrname(rec.tid)]
            return ('E', int(cL - rec.pos))
        else: 
            # Read is forward wrt contig
            return ('B', int(rec.positions[-1])+1)
    ##############################################################################


    r1Func = readSense if r1IsForward else readAntisense
    r2Func = readAntisense if r2IsForward else readSense

    for r1, r2 in yieldPairsSmart(bamFile):
        pairCount += 1

        if r1.tid == r2.tid:

        # Case 1: The reads align to the same contig
            r1F =  not r1.is_reverse
            r2F =  not r2.is_reverse

            # We must check that the read pairs aligning to the same contig do so concordantly.
            # If the alignment orientation of read 1 (or 2) matches the expectation, then read 1 (or 2) says contig is forward. 
            # So orientationOK determines if read1 and read2 agree on the orientation of the contig.
            # If the contig is forward (w.r.t insert), we expect read 1 to the left of read 2 in the contig (and vice versa is contig is reverse).
            # So placementOK determines if read1 and read2 have the correct placement within the contig

            orientationOK = ((r1F == r1IsForward) == (r2F == r2IsForward)) # Check on relative orientation
            placementOK = ((r1F == r1IsForward) == (r1.pos < r2.pos) == (r1.positions[-1] < r2.positions[-1])) # Check on relative placement
            if  (orientationOK and placementOK):
                H.add(abs(r1.tlen))
                concordantCount += 1
                continue

        # Case 2: The read pairs align to different contigs, or
        # they align discordantly to the same contig
        # Create a link between contigs
        or1, d1 = r1Func(r1)
        or2, d2 = r2Func(r2)
        assert(d1>0 and d2>0)
        contig1 = bamFile.getrname(r1.tid)
        contig2 = bamFile.getrname(r2.tid)

        if contig1 < contig2:
            linkDict[(contig1, contig2)].append((or1, d1, or2, d2))
        elif contig2 < contig1:
            linkDict[(contig2, contig1)].append((or2, d2, or1,d1))
        elif contig1 == contig2:
            evidence = (or1, d1, or2, d2) if or1 < or2  else (or2, d2, or1, d1)
            linkDict[(contig1, contig2)].append(evidence)
        else:
            assert(False)
        linkCount += 1

    print 'Read %i pairs. Made %i links'%(pairCount, linkCount)
    print 'Insert Size Distribution:\nMean: %10.2f\nMedian: %10.2f\nstd: %10.2f'%(H.mean(), H.median(), H.std())

    # Pack a dictionary object with results to return
    res = {}
    res['histogram'] = H
    res['insMean'] = H.mean()
    res['insMedian'] = H.median()
    res['insStd'] = H.std()
    res['mappedPairs'] = pairCount
    res['concordantPairs'] = concordantCount
    res['numLinks'] = linkCount
    res['linkDict'] = linkDict
    return res

###############################################
# Take a linkDict generated by readBam_ver2
# and bundle links between contigs to produce a 
# distance estimate of the gap between contigs.
# A negative gap implies an overlap
#
# minStd is the minimum standard deviation employed for gap estimate (in bp)
def makeBundles(linkDict, insMean, insStd, maxDev=5, minStd=5.0):
    maxOverlap = insMean + maxDev*insStd

    bundle_fields = ['id', 'n1','n2','e1','e2','m','std','n', 'd1', 'd2']
    # id: unique id for the bundle
    # n1: node 1
    # n2: node 2
    # e1: end of node 1 
    # e2: end of node 2
    # m: gap estimate mean
    # std: gap estiamte std
    # n: number of pairs (links) used to create bundle
    # d1: number of insert bases from read 1 to end of node 1
    # d2: number of insert bases from read 2 to end of node 2
    bundleDict = defaultdict(list)                                                                                                                                                           
    bid = 0
    badOverlapCount = 0 # Number of links that imply an overlap which is too large.
    for (n1,n2), evidenceList in linkDict.iteritems():

        # Group the links based on the relative orientation they imply between n1 & n2
        linkGroups = defaultdict(list)
        for (e1,d1,e2,d2) in evidenceList:
            gap = insMean - (d1+d2)
            if (gap >= 0) or (abs(gap) < maxOverlap):
                linkGroups[(e1,e2)].append((gap,d1,d2))
            else:
                badOverlapCount += 1

        # For each relative orientation, make an estimate of the gap
        # Note: This estimate could be made more robust!
        for (e1,e2), links in linkGroups.iteritems():
            gaps, d1, d2 = zip(*links)
            gaps = np.array(gaps)
            n = gaps.shape[0]
            m = np.mean(gaps)
            std = max(insStd / np.sqrt(n), minStd)
            b = dict(zip(bundle_fields, [bid, n1,n2,e1,e2,m,std,n,d1,d2]))
            bundleDict[(n1,n2)].append(b)
            bid += 1

    # Pack a dictionary with the results
    res = {}
    res['bundleDict'] = bundleDict
    res['badLinks'] = badOverlapCount
    res['numBundles'] = bid
    return res

    
#################################################
# Summarize the links and bundle information
# resLinks: result from readBam_ver2
# resBundles: result returned from makeBundles
def writeSummary(resLinks, resBundles, foutName):
    fout = open(foutName, 'w')
    sList = []
    sList.append('Mapped Pairs: %i'%resLinks['mappedPairs'])
    p = 100*float(resLinks['concordantPairs'])/resLinks['mappedPairs']
    sList.append('Concordant Pairs: %i (%6.3f %%)'%(resLinks['concordantPairs'], p))
    sList.append('Number of Links: %i'%resLinks['numLinks'])
    sList.append('Insert mean: %f'%resLinks['insMean'])
    sList.append('Insert median: %f'%resLinks['insMedian'])
    sList.append('Insert std.: %f'%resLinks['insStd'])
    sList.append('Number of Bundles: %i'%resBundles['numBundles'])
    p = 100*float(resBundles['badLinks'])/resLinks['numLinks']
    sList.append('Number of Bad Links (overlap too large): %i (%6.3f %%)'%(resBundles['badLinks'],p))
    fout.write('\n'.join(sList))
    fout.close()

#################################################
# Convert a bundle to a string
def bundle2Str(b):
    kList = ['n1', 'e1', 'n2', 'e2', 'm', 'std', 'n']
    d1S = ','.join(map(str,b['d1']))
    d2S = ','.join(map(str,b['d2']))
    strFields = []
    strFields.append('bundle-%i'%b['id'])
    strFields.extend(str(b[k]) for k in kList)
    strFields.append(d1S)
    strFields.append(d2S)
    outS = ' '.join(f for f in strFields)
    return outS

#################################################
# Convert a string to a bundle
def str2bundle(s):
    fields = s.strip().split()
    d1 = [int(v) for v in fields[-2].split(',')]
    d2 = [int(v) for v in fields[-1].split(',')]
    fields = fields[:-2]
    fieldNames = ['id', 'n1', 'e1', 'n2', 'e2', 'm', 'std', 'n']
    fieldTypes = [str, str, str, str, str, float, float, int]
    b = dict((fn, ft(fv)) for fn,ft,fv in zip(fieldNames, fieldTypes, fields))
    b['d1'] = d1
    b['d2'] = d2
    return b

#################################################
def writeBundles(bundleDict, outFileName):
    fout = open(outFileName, 'w')
    # bundleDict keys are (contig1, contig) pairs, and values are list of bundles between them
    # Here, we flatten this into a new dictionary based on bundle id
    bd = dict((b['id'], b) for bList in bundleDict.itervalues() for b in bList)
    kSorted = sorted(bd.keys())
    for k in kSorted:
        fout.write('%s\n'%bundle2Str(bd[k]))
    fout.close()

###############################################
# Read SamFile corresponding to alignment to contig ends
# This function expects contig end to have a very specific name, which
# encodes the contig end 
# Determine how many pairs are in same contig.
# Make histogram for insert size (controlled by insLow, insHigh, insInc)
# Generate links between different contigs
def readSamEnds(samFileName, fileType='sam'):

    def decodeContigName(contigName):
        fields = contigName.split('.')
        contigId,coords,orientation = fields
        coords = tuple(map(int,coords.split(':')))
        assert(len(coords)==2)
        return (contigId,coords,orientation)


    if fileType == 'bam':
        samFile = pysam.Samfile(samFileName, 'rb')
    else:
        samFile = pysam.Samfile(samFileName, 'r')

    contigLengthDict = dict(zip(samFile.references, samFile.lengths))

    invalidCount = 0
    pairCount = 0
    diffContigs = 0
    sameContigs = 0
    links = [] # Links between different contigs
    for r1, r2 in yieldPairsSmart(samFile):

        pairCount += 1
        if r1.tid == r2.tid:
            sameContigs += 1
            continue

        diffContigs += 1
        # Create a link between contigs
        # Assume read1 is forward and read2 is reverse
        # Therefore, if read1 is mapped as forward to contig1, then contig1 is forward
        # and if read2 is mapped as forward to contig2, then contig2 is reverse.
        contig1 = samFile.getrname(r1.tid)
        contig1L = contigLengthDict[contig1]
        contig1Name, contig1Coords, contig1Or = decodeContigName(contig1)
        contig2 = samFile.getrname(r2.tid)
        contig2L = contigLengthDict[contig2]
        contig2Name, contig2Coords, contig2Or = decodeContigName(contig2)
        contig2 = samFile.getrname(r2.tid)

        # Check for proper orientation of reads
        r1Invalid = ((r1.is_reverse) and (contig1Or=='R')) or ((not r1.is_reverse) and (contig1Or=='L'))
        r2Invalid = ((r2.is_reverse) and (contig2Or=='R')) or ((not r2.is_reverse) and (contig2Or=='L'))
        if r1Invalid or r2Invalid:
            invalidCount += 1
            continue

        # d1 and d2 capture the amount of the insert size that are contained in each contig.
        # Ex:
        # |-------------------|
        #       r1 |--->
        #          |<----d1-->|

        # We assume that the orientation of mates is FR:

        # |----> read 1
        # |-------------------->| template
        # |<--------------------| template
        #         read 2   <----| 

        if not r1.is_reverse:
            # Read 1 is forward with respect to its contig
            # This implies link to end of contig
            # |------------------> contig
            #               |-->  read
            end1 = 'E'
            d1 = int(contig1L - r1.pos) 
        else:
            # Read 1 is reverse with respect to its contig
            # This implies link to beginning of contig
            # |------------------> contig
            #   <---|  read
            end1 = 'B'
            d1 = int(r1.positions[-1])+1

        if not r2.is_reverse:
            # Read 1 is forward with respect to its contig
            # This implies link to end of contig
            # |------------------> contig
            #               |-->  read
            end2 = 'E'
            d2 = int(contig2L - r2.pos)
        else:
            end2 = 'B'
            d2 = int(r2.positions[-1])+1
            # Read 1 is reverse with respect to its contig
            # This implies link to beginning of contig
            # |------------------> contig
            #   <---|  read
        links.append((contig1Name, end1, d1, contig2Name, end2, d2))

    print 'Read %i pairs.\n%i pairs are invalid.\n%i pairs are in different contigs.'%(pairCount, invalidCount, diffContigs)
    return links

###############################################
# Create a dictionary where linkDict[node1][node2] gives a list of links between node1 & node2
# The output link is given as a tuple : (node 1 end, node 2 end, separation)
# where separation is the MLE for separation between the ends of nodes for that single link.
# Inputs:
# Links: List of (contig1, or1, d1, contig2, or2, d2) (computed by readSamEnds)
# Sep is the mean separation between read pairs for this list of links.
def gatherLinks(sep, links):
    linkDictMaker = lambda : defaultdict(list)
    linkDict = defaultdict(linkDictMaker)
    #links.append((contig1, or1, d1, contig2, or2, d2))
    reverseOr = lambda s: 'F' if s == 'R' else 'R'
    for l in links:
        contig1, end1, d1, contig2, end2, d2 = l
        linkSep = sep - (d1+d2)
        linkDict[contig1][contig2].append((end1, end2, linkSep))
        linkDict[contig2][contig1].append((end2, end1, linkSep))
    return linkDict

###############################################
def contigCoverage(samFileName, insMean, insStd, fileType='sam'):
    contigReads = defaultdict(list) #Key: ContigId Val: List of tuples (start, end) for read locations)
    contigInserts = defaultdict(list) #Key: ContigId Val: List of tuples (start, end) for insert locations
    if fileType == 'bam':
        samFile = pysam.Samfile(samFileName, 'rb')
    else:
        samFile = pysam.Samfile(samFileName, 'r')

    contigLengthDict = dict(zip(samFile.references, samFile.lengths))
    # Sort contigs by length
    sortedContigs = sorted((length,contigId) for contigId,length in contigLengthDict.iteritems())[::-1]

    # Read paired reads from SamFile
    sys.stdout.write('Processing file %s....'%samFileName)
    sys.stdout.flush()
    for r1,r2 in yieldPairs(samFile):
        assert(r1.is_reverse != r2.is_reverse)
        read1 = r1 if not r1.is_reverse else r2
        read2 = r1 if r1.is_reverse else r2
        assert(read1.rname == read2.rname)
        contigId = read1.rname
        contigReads[contigId].append((int(read1.positions[0]), int(read1.positions[-1])))
        contigReads[contigId].append((int(read2.positions[0]), int(read2.positions[-1])))
        contigInserts[contigId].append((int(read1.positions[0]), int(read2.positions[-1])))
    sys.stdout.write('DONE!\n') 
    sys.stdout.flush()
    # TO DO: for each contig, turn contigReads and contigInserts into numpy array. Sort. Then for each position in the contig,
    # produce a value for read coverage, insert coverage, and mean insert sizes.
    assert(sorted(contigReads.keys()) == sorted(contigInserts.keys()))
    fout = open('contigData.txt', 'w')
    numContigs = len(contigReads)
    contigCount = 0
    for contigL, contigId in sortedContigs:
        tid = samFile.gettid(contigId)
        sys.stdout.write('Processing Contig %i of %i (%f percent)\n'%(contigCount, numContigs,100.0*contigCount/numContigs))
        sys.stdout.flush()
        contigCount += 1
        fout.write('-'*50 + '\n')
        fout.write('Contig: %s\n'%contigId)
        #contigId = samFile.getrname(tid)
        #contigL = contigLengthDict[contigId]

        # Get read locations (start & end pts) aligned to this contig, and sort
        readLocs = np.array(contigReads[tid])
        numReads = readLocs.shape[0]
        readSortInd = np.argsort(readLocs, 0)

        # Get insert locations (start & end pts) aligned to this contig and sort
        insertLocs = np.array(contigInserts[tid])
        numInserts = insertLocs.shape[0]
        insertSortInd = np.argsort(insertLocs, 0)

        assert (numReads == 2 * numInserts)
   
        # Pack the start and end locations into an array, with the start or end loc and it's sorted rank
        readStart = np.zeros((numReads,2))
        readEnd = np.zeros((numReads, 2))
        if numReads > 0:
            readStart[:,0] = readLocs[readSortInd[:,0],0]
            readStart[:,1] = readSortInd[:,0]
            readEnd[:,0] = readLocs[readSortInd[:,1],1]
            readEnd[:,1] = readSortInd[:,1]

            insertStart = np.zeros((numInserts,2))
            insertStart[:,0] = insertLocs[insertSortInd[:,0],0]
            insertStart[:,1] = insertSortInd[:,0]
            insertEnd = np.zeros((numInserts, 2))
            insertEnd[:,0] = insertLocs[insertSortInd[:,1],1]
            insertEnd[:,1] = insertSortInd[:,1]


        # Annotate every position in contig by it's read coverage, insert coverage, and C/E statistic
        # Write these to file
        curFragments = {}
        curReadCov = 0
        curInsCov = 0
        rsInd = 0
        rs = readStart[rsInd, 0] if rsInd < numReads else contigL# The next read start position
        reInd = 0
        re = readEnd[reInd, 0] if reInd < numReads else contigL # The next read end position
        insSInd = 0
        insS = insertStart[insSInd, 0] if insSInd < numReads else contigL# The next insert start position
        insEInd = 0
        insE = insertEnd[insEInd,0] if insEInd < numReads else contigL # The next insert end position
        curCE = 0.0
        for curLoc in range(contigL):
            # check new reads that start at this location
            assert(curLoc <= rs)
            while curLoc == rs:
                curReadCov += 1
                rsInd += 1
                if rsInd < numReads:
                    assert readStart[rsInd, 0] >= rs
                rs = readStart[rsInd, 0] if rsInd < numReads else contigL
            # check reads that end at this location
            assert (curLoc <= re + 1)
            while curLoc == re + 1:
                curReadCov -= 1
                reInd += 1
                if reInd < numReads:
                    assert readEnd[reInd, 0] >= re
                re = readEnd[reInd, 0] if reInd < numReads else contigL
            # new inserts that start at this location
            assert (curLoc <= insS)
            recomputeCE = False
            while curLoc == insS:
                recomputeCE = True
                unsortedInd = insertStart[insSInd,1]
                insSize = insertLocs[unsortedInd,1] - insertLocs[unsortedInd,0] + 1
                assert (unsortedInd not in curFragments)
                curFragments[unsortedInd] = insSize
                curInsCov += 1
                insSInd += 1
                if insSInd < numInserts:
                    assert insertStart[insSInd, 0] >= insS
                insS = insertStart[insSInd, 0] if insSInd < numInserts else contigL

            # inserts that end at this location
            while curLoc == insE + 1:
                recomputeCE = True
                unsortedInd = insertEnd[insEInd,1]
                curFragments.pop(unsortedInd)
                curInsCov -= 1
                insEInd += 1
                if insEInd < numInserts:
                    assert insertEnd[insEInd, 0] >= insE
                insE = insertEnd[insEInd,0] if insEInd < numInserts else contigL
            if recomputeCE:
                meanInsSpan = np.mean(curFragments.values())
                curCE = (meanInsSpan - insMean) / (insStd / np.sqrt(curInsCov)) if curInsCov > 0 else 0.0
            fout.write('%i\t%i\t%i\t%f\n'%(curLoc, curReadCov, curInsCov, curCE))
    fout.close() 
