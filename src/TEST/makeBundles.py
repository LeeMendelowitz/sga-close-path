import investigateSam as invS

bamFile = 'assemble.K27.X2.m70-contigs.alignments.bam'
#bamFile = 'assemble.K27.X2.m70-contigs.alignments.namesort.head.sam'
insLow = 0
insHigh = 500
insInc = 5
print 'Creating links from bam file %s...'%bamFile
#resLinks = invS.readBam_ver2(bamFile, insLow, insHigh, insInc, fileType='sam')
resLinks = invS.readBam_ver2(bamFile, insLow, insHigh, insInc, fileType='bam')
linkDict = resLinks['linkDict']

print 'Creating bundles from links...'
resBundles = invS.makeBundles(linkDict, resLinks['insMean'], resLinks['insStd'])
bundleDict = resBundles['bundleDict']

print 'Writing bundle file...'
bundleFileOut = 'assemble.K27.X2.m70-contigs.alignments.bundles3'
summaryFileOut = 'assemble.K27.X2.m70-contigs.alignments.summary'
invS.writeBundles(bundleDict, bundleFileOut)

invS.writeSummary(resLinks, resBundles, summaryFileOut)


