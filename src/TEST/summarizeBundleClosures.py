# Parse the output of the SGA Path Closer program and summarize
# the bundles.

# TO DO: USe pandas to parse?

def readStatusFile(fname):
    fin = open(fname)
    fin.next() # skip header
    #statusDict = {}
    statusRecs = []

    str2bool = lambda v: not (v=='0')
    fieldNames = ['id', 'numClosures', 'tooRepetative', 'overlapTooLarge']
    fieldTypes = [str, int, str2bool, str2bool]

    for l in fin:
        fields = l.strip('\n').split('\t')
        bDict = dict((fn, ft(fv)) for fn,ft,fv in zip(fieldNames, fieldTypes, fields))
        #statusDict[bDict['id']] = bDict
        statusRecs.append(bDict)

    fin.close()

    return statusRecs

def readStatsFile(fname):
    fin = open(fname)
    fin.next() # skip header
    statsRecs = []

    st2bool = lambda v: not (v=='0')
    def strToList(s):
            
            return [int(v) for v in s.split(',') if v]

    fieldNames = ['id', 'gapEst', 'stdEst', 'numLinks', 'numClosures', 'closureLengths']
    fieldTypes = [str, float, float, int, int, strToList]

    for l in fin:
        fields = l.strip('\n').split('\t')
        bDict = dict((fn, ft(fv)) for fn,ft,fv in zip(fieldNames, fieldTypes, fields))
        statsRecs.append(bDict)
        #statsDict[bDict['id']] = bDict
    fin.close()
    return statsRecs

def makeSummaryDict(statsFile, statusFile):
    statusRecs = readStatusFile(statusFile)
    statsRecs = readStatsFile(statsFile)
    assert(len(statusRecs)==len(statsRecs))
    summaryDict = {}
    for status, stats in zip(statusRecs, statsRecs):
        assert(status['id']==stats['id'])
        summaryDict[status['id']] = dict(status.items() + stats.items())
    return summaryDict

def countLinks(summaryDict):
    linksNoPath = 0
    linksUniquePath = 0
    linksMultiPath = 0
    for s, bd in summaryDict.iteritems():
        if bd['numClosures'] == 0:
            linksNoPath += bd['numLinks']
        elif bd['numClosures'] == 1:
            linksUniquePath += bd['numLinks']
        elif bd['numClosures'] > 1:
            linksMultiPath += bd['numLinks']
        else:
            assert(False)
    return (linksNoPath, linksUniquePath, linksMultiPath)

    

if __name__ == '__main__':
    test()
