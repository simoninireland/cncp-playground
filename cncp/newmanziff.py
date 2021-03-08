from epyc import Lab, Experiment, ResultsDict
from epydemic import NetworkExperiment, Node, Edge
from networkx import Graph
import numpy

from typing import Any, Dict, Union, Iterable, List, Final, Optional, cast


class NewmanZiff(NetworkExperiment):
    '''Base class for the Newman-Ziff site and bond percolation algorithm.
    
    :param g: (optional) the underlying network or generator
    :param samples: (optional) number of samples or list of sample points (defaults to 100)'''

    def __init__(self, g : Graph =None, samples : Union[int, Iterable[float]] =None):
        super(NewmanZiff, self).__init__(g)

        # fill in default
        if samples is None:
            samples = 100
        if isinstance(samples, int):
            samples = numpy.linspace(0.0, 1.0, num=samples, endpoint=True)
        self._samples : List[float] = sorted(numpy.unique(list(cast(Iterable[float], samples))))

        # components data structure is initially empty
        self._components : numpy.ndarray = None

        # no component
        self._gcc = 0

    def tearDown(self):
        '''Throw away the components data structure at tear-down.'''
        self._components = None
        super(NewmanZiff, self).tearDown()

    def report(self, params : Dict[str, Any], meta : Dict[str, Any], res : Union[Dict[str, Any], List[ResultsDict]]) -> ResultsDict:
        for rc in res:
            rc[Experiment.PARAMETERS] = params.copy()
            rc[Experiment.METADATA] = meta.copy()
        return res

    
class BondPercolation(NewmanZiff):
    '''A bond (edge) percolation experiment. This experiment computes the size of
    the giant connected component (GCC) as edges are "occupied" within the underlying network.
    It samples the GCC at a given sequence of occupation probabilities, returning
    a time series of the growth of the GCC.

    :param g: (optional) the underlying network or generator
    :param samples: (optional) number of samples or list of sample points (defaults to 100)'''

    # Synthesised parameters
    P : Final[str] = 'epydemic.bondpercolation.pOccupied'    #: Parameter holding percolation threshold.

    # Experimental results
    GCC : Final[str] = 'epydemic.bondpercolation.gcc'        #: Result holding size of GCC.

    def __init__(self, g : Graph =None, samples : Union[int, Iterable[float]] =None):
        super(BondPercolation, self).__init__(g, samples)

    def setUp(self, params : Dict[str, Any]):
        '''Set up the process, creating the initial components data structure from the
        underlying network.

        :param params: the experimental parameters'''
        super(BondPercolation, self).setUp(params)
        self._components = numpy.full(N, -1, numpy.int32)
        self._gcc = 1   # initially all nodes are individual components, unconnected by occupied edges

    def rootOf(self, n : Node) -> Node:
        np = self._components[n]
        if np < 0:
            # n is the root, return it
            return n
        else:
            # n has a parent, follow the tree to it
            r = self.rootOf(np)

            # update our component record to point to the root
            self._components[n] = r

            # return the root
            return r    

    def occupy(self, n : Node, m : Node) -> Optional[int]:
        nr = self.rootOf(n)
        mr = self.rootOf(m)
        if mr != nr:
            # nodes are in different components, join them together
            return self.join(nr, mr)
        else:
            # nodes are in the same component, do nothing
            return None
        
    def join(self, c1 : Node, c2 : Node) -> int:
        # extract the size of the second compooent
        msize = self._components[c2]
        
        # join the second compoent to the first 
        self._components[c2] = c1
        
        # update the size of the first component
        self._components[c1] += msize
        
        # update the GCC
        self._gcc = max(self._gcc, -self._components[c1] )
        
        # return the size of the new component
        return -self._components[c1]        
    
    def sample(self, p : float) -> Dict[str, Any]:
        res = dict()
        res[self.P] = p        
        res[self.GCC] = self._gcc        
        return res
    
    def percolate(self, es : List[Edge]) -> List[Dict[str, Any]]:
        # take an initial sample if requested
        samples = []
        samplePoint = 0
        if self._samples[samplePoint] == 0.0:
            samples.append(self.sample(self._samples[samplePoint]))
            samplePoint += 1

        # percolate the network
        M = len(es)
        for i in range(M):
            (n, m) = es[i]

            # occupy the edge
            csize = self.occupy(n, m)

            # take a sample if this is a sample point
            if  (i + 1) / M >= self._samples[samplePoint]:
                # we're at the closest probability after the requested sample point,
                # so build the sample
                samples.append(self.sample(self._samples[samplePoint]))

                # if we've collected all the samples we want, bail out
                samplePoint += 1
                if samplePoint > len(self._samples):
                    break
                    
        return samples
    
    def do(self, params : Dict[str, Any]) -> List[ResultsDict]:
        '''Perform the bond percolation process.

        :param params: experimental parameters
        :returns: a list of dicts of experimental results'''        
        # extract and shuffle the edges
        g = self.network()
        es = list(g.edges()).copy()
        numpy.random.shuffle(es)

        res = self.percolate(es)
        
        rcs = []
        for obs in res:
            rc = self.resultsdict()
            rc[Experiment.RESULTS] = obs
            rcs.append(rc)
        return rcs
    
