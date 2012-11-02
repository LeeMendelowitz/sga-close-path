import numpy as np

###############################################################################
# Histogram class for storing counts
class Histogram(object):

    def __init__(self, binEdges):
        '''
        Generate a Histogram instance with bin boundaries in bins.
        The intervals are closed on the left. Ex: bins = [0, 5, 10]
        Gives three intervals: [0,5) and [5,10) and [10,infinity)
        '''
        self.binEdges = np.sort(binEdges)
        self.numBins = self.binEdges.shape[0]+1 # Include the boundary bins for data that is too small or too large
        self.rBinEdges = self.binEdges[1:]
        self.counts = np.zeros(self.numBins)

    def add(self, nums):
        '''
        Add a float (or a list of floats) into the histogram
        '''
        if (isinstance(nums, float) or isinstance(nums, int)):
            nums = [nums]
        # Do binary search to find the right bin
        inds = np.searchsorted(self.binEdges, nums, side='right')
        for ind in inds:
            self.counts[ind] += 1


    def __str__(self):
        outS = ''
        outS += '<%f: %i\n'%(self.binEdges[0], self.counts[0])
        outS += ''.join('[%f, %f): %i\n'%(self.binEdges[i-1], self.binEdges[i], self.counts[i]) for i in range(1,self.numBins-1))
        outS += '>=%f: %i\n'%(self.binEdges[-1], self.counts[-1])
        return outS

    def mean(self):
        n = np.sum(self.counts[1:-1])
        binCenters = 0.5*(self.binEdges[0:-1] + self.binEdges[1:])
        p = self.counts[1:-1]/float(n)
        return np.sum(p*binCenters)

    def var(self):
        binCenters = 0.5*(self.binEdges[0:-1] + self.binEdges[1:])
        n = np.sum(self.counts[1:-1])
        p = self.counts[1:-1]/float(n)
        var = np.sum(p*(binCenters**2)) - (np.sum(p*binCenters))**2
        return var

    def std(self):
        return np.sqrt(self.var())

    def median(self):
        sums = np.cumsum(self.counts[1:-1])
        binCenters = 0.5*(self.binEdges[0:-1] + self.binEdges[1:])
        assert ( sums.shape == binCenters.shape)
        n = int(sums[-1])
        n_over_2 = n/2.0
        if n%2 == 0:
            ind = np.nonzero(sums >= n_over_2)[0][0]
            return binCenters[ind]
        else:
            n_over_2_left = n_over_2 - 0.5
            n_over_2_right = n_over_2 + 0.5
            indLeft = np.nonzero(sums >= n_over_2_left)[0][0]
            indRight = np.nonzero(sums >= n_over_2_right)[0][0]
            return 0.5*(binCenters[indLeft] + binCenters[indRight])
