from epyc import Lab, Experiment, ResultsDict
from epydemic import NetworkExperiment, Node, Edge
from networkx import Graph
import numpy
from cncp import NewmanZiff

from typing import Any, Dict, Union, Iterable, List, Final, Optional, cast

class ResidualBondPercolation(NewmanZiff):

    # Experimental results
    P_RESIDUAL_STEM : Final[str] = 'epydemic.residualbondpercolation.pOccupied'
    DEPTH_RESIDUAL : Final[str] = 'epydemic.residualbondpercolation.depth'
    N_RESIDUAL : Final[str] = 'epydemic.residualbondpercolation.N'
    M_RESIDUAL : Final[str] = 'epydemic.residualbondpercolation.M'
    GCC_RESIDUAL : Final[str] = 'epydemic.residualbondpercolation.gcc'
        
    @classmethod
    def P_RESIDUAL(cls, depth : int) -> str:
        '''Return the experimental result corresponding to the occupation probability
        at the given depth: 0 for the primary network, 1 for the first residual network,
        and so on.
        
        :param residual: the network depth
        :returns: the experimental result tag'''
        return '{stem}-{l}'.format(stem=cls.P_RESIDUAL_STEM, l=depth)
    
    def __init__(self, g : Graph =None, samples : Union[int, Iterable[float]] =None, residuals =1):
        super(ResidualBondPercolation, self).__init__(g, samples)
        self._residuals = residuals
    
    def setUp(self, params : Dict[str, Any]):
        super(ResidualBondPercolation, self).setUp(params)
        self._networkIndex = 0
        self._parent = 0
        self._phis = dict()
        self._gcc = 1   # initially all nodes are individual components, unconnected by occupied edges
        
        # all nodes are singleton components in the primary network
        N = self.network().order()
        self._components = numpy.full(N, -1, numpy.int32)
        self._networks = numpy.full(N, 1, numpy.int16)
        
    def rootOf(self, n : Node, network : int) -> Optional[Node]:
        np = self._components[n]
        nn = self._networks[n]
        
        # for a node to be available to this percolation process it has
        # to be one of:
        # - an unoccupied element of the base network 0
        # - an occupied element of the network we're currently building
        # - an element of a network that's later than out parent network
        # These last two conditions overlap (we're later than our parent too),
        # but in the third case the node is treated as an unoccupied singleton
        # and not as a member of a component.
        
        if np == -1:
            # node is an unoccupied singleton, move it to our network and return it
            self._networks[n] = network
            return n
        elif nn == network:
            # node is in our network, use it
            if np < 0:
                # node is a component root, return it
                return n
            else:
                # node is an interior node of a component, track it to the root
                r = self.rootOf(np, network)
                #assert(r is not None)
                #assert(self._networks[r] == self._networks[n])
                
                # short-cut to the root
                self._components[n] = r
                
                # return the root
                return r
        elif nn > self._parent:
            # node is in another available network, reset it as a singleton in ours and return
            self._components[n] = -1
            self._networks[n] = network
            return n
        else:
            # node is in a network inaccessible to us, ignore
            return None

    def orderOfNetwork(self, network : int) -> int:
        if network < 1:
            return 0
        else:
            return numpy.count_nonzero(self._networks == network)

    def occupy(self, n : Node, m : Node, network) -> Optional[int]:
        nr = self.rootOf(n, network)
        mr = self.rootOf(m, network)
        if nr is None or mr is None:
            # not a node in the right network, ignore
            return None
        elif mr != nr:
            # nodes are in different components in this network, join them together
            #assert(self._networks[nr] == self._networks[mr])
            return self.join(nr, mr)
        else:
            # nodes are in the same component, do nothing
            #assert(self._networks[nr] == self._networks[mr])
            return None
        
    def join(self, c1 : Node, c2 : Node) -> int:
        # extract the sizes of the compoents
        c1size = -self._components[c1]
        c2size = -self._components[c2]
        
        # join the second component to the first 
        self._components[c2] = c1
        
        # update the size of the first component
        self._components[c1] -= c2size
        
        # update the running GCC size
        csize = c1size + c2size
        self._gcc = max(self._gcc, csize)
        
        # return the size of the new component
        return csize

    def sample(self, p, es, nexti, N, M, depth, network):
        samples = []
        
        # take the sample at this point
        phi = ResidualBondPercolation.P_RESIDUAL(depth)
        res = self._phis.copy()
        res[self.DEPTH_RESIDUAL] = depth
        res[self.N_RESIDUAL] = N
        res[self.M_RESIDUAL] = M
        res[phi] = p
        res[self.GCC_RESIDUAL] = self._gcc
        samples.append(res)
        
        # if we're re-percolating, do that to get more samples
        if depth < self._residuals:
            es_residual = es[nexti:]
            samples.extend(self.repercolate(p, es_residual, depth, network))
        
        return samples

    def repercolate(self, p, es, depth, network):
        phi = ResidualBondPercolation.P_RESIDUAL(depth)
        self._phis[phi] = p
        oldparent = self._parent
        oldgcc = self._gcc
        self._parent = network
        self._gcc = 1
        ss = self.percolate(es, depth + 1)
        del self._phis[phi]
        self._gcc = oldgcc
        self._parent = oldparent        
        return ss
        
    def percolate(self, es, depth):
        N = self.network().order() - self.orderOfNetwork(self._parent)
        M = len(es)
        self._networkIndex += 1
        network = self._networkIndex
        samples = []
        samplePoint = 0

        # take an initial sample if requested
        if self._samples[samplePoint] == 0.0:
            samples.extend(self.sample(self._samples[samplePoint], es, 0, N, M, depth, network))
            samplePoint += 1
            
        # percolate the network
        for i in range(M):
            # if we've collected all the samples we want, bail out
            if samplePoint >= len(self._samples):
                break
                    
            # occupy the edge
            (n, m) = es[i]
            csize = self.occupy(n, m, network)

            # take a sample if this is a sample point
            if  (i + 1) / M >= self._samples[samplePoint]:
                # we're at the closest probability after the requested sample point,
                # so build the sample
                samples.extend(self.sample(self._samples[samplePoint], es, i + 1, N, M, depth, network))
                samplePoint += 1
                    
        return samples

    def do(self, params : Dict[str, Any]) -> List[Dict[str, Any]]:
        '''Perform the bond percolation process.

        :param params: experimental parameters
        :returns: a list of dicts of experimental results'''      
        # extract and shuffle the edges
        g = self.network()
        es = list(g.edges()).copy()
        numpy.random.shuffle(es)
        
        res = self.percolate(es, 0)
    
        rcs = []
        for obs in res:
            rc = self.resultsdict()
            rc[Experiment.RESULTS] = obs
            rcs.append(rc)
        return rcs
    
