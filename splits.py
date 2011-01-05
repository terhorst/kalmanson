# Parallel code
from sage.all import Set, combinations_iterator, binomial
import operator
import itertools as it
from memoized import memoized
from kalmanson import ray_to_splits
from pq_split_system import PQSplitSystem as SplitSystem

class Split(object):
    def __init__(self, X, A):
        A = Set(A)
        X = Set(X)
        assert A <= X
        self.X = X
        self.A = A if min(X) in A else X - A
        self.B = X - self.A

    def __cmp__(self, other_split):
        return self.blocks().__cmp__(other_split.blocks())

    def __hash__(self):
        return self.blocks().__hash__()

    def __getstate__(self):
        return {'X': self.X, 'A': self.A, 'B': self.B}

    def __setstate__(self, state):
        self.X = state['X']
        self.A = state['A']
        self.B = state['B']

    def __repr__(self):
        return str(self.blocks())

    def merge(self, other):
        return Split(self.X, self.A.intersection(other.A))

    def blocks(self):
        return Set([self.A, self.X - self.A])

def all_splits(n):
    X = Set(range(1, n+1))
    return Set(Split(X,b) for i in range(2, X.cardinality() - 1) \
            for b in combinations_iterator(X, i))

def check_circular(ss):
    return [ss, ss.is_circular()]

@memoized
def circular_splits(n, k):
    from sage.parallel.ncpus import ncpus
    from multiprocessing import Pool
    from sage.misc.fpickle import pickle_function, call_pickled_function
    from twisted.internet import reactor   # do not delete this (!)  -- see trac 8785

    allsp = all_splits(n)

    if k==1:
        return [SplitSystem([s]) for s in allsp]

    Snk1 = circular_splits(n, k-1)
    Snk = (ss.add_split(sp) for ss in Snk1 for sp in allsp if not ss.contains(sp))

    p = Pool(ncpus())
    chunksize = min(1000, int(len(Snk1) * len(allsp) / ncpus()))
    res = p.imap_unordered(check_circular, Snk, chunksize)
    return Set(ss for ss,circ in res if circ)

def fvector(n):
    return [len(circular_splits(n,k)) for k in range(1, n*(n-3)/2 + 1)]

def makesplit(n, lst):
    lst = lst if 1 in lst else Set(range(1, n+1)) - Set(lst)
    return reduce(operator.or_, [2**(n - i) for i in lst])

def split_system_from_cone(n, cone):
    X = range(1,n+1)
    return SplitSystem(Split(X, ray_to_splits(n,r)) for r in cone.rays())

def splits2hs(n, ss):
    ret = []
    tpl = "fromList [%s]"
    for sps in ss:
        ret.append(tpl % ",".join([str(makesplit(n, sp.A)) for sp in sps.splits()]))
    return tpl % ", ".join(ret)
