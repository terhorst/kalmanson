from sage.all import Set, combinations_iterator, parallel, binomial
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


@parallel(p_iter='multiprocessing')
def __circular_splits_helper(sp):
    return sp.is_circular()

@parallel(p_iter='multiprocessing')
def __weak_splits_helper(sp):
    from splits import SplitSystem
    return SplitSystem(sp).is_weakly_compatible()

def weakly_compatible_splits(n, k):
    return __splits_checker(n, k, __weak_splits_helper)

def circular_splits(n, k):
    return __splits_checker(n, k, __circular_splits_helper)

def fvector(n):
    return [len(circular_splits(n,k)) for k in range(1, n*(n-3)/2 + 1)]

@memoized
def __splits_checker(n, k, helper):
    allsp = all_splits(n)
    S = [SplitSystem([s]) for s in allsp]
    if k==1:
        return S
    Snm1 = __splits_checker(n, k-1, helper)
    other_ss = Set(SplitSystem(sp) for sp in combinations_iterator(allsp, k-1)) - Snm1
    impossible_ss = Set(ss.add_split(sp) for ss,sp in it.product(other_ss, allsp))
    if binomial(n, k) > ( len(S) * len(Snm1) ):
        ss = filter(lambda ss: len(ss)==k and ss not in impossible_ss,
                (ss1.join(ss2) for ss1,ss2 in it.product(S, Snm1)))
    else:
        ss = filter(lambda ss: ss not in impossible_ss,
                (SplitSystem(sp) for sp in combinations_iterator(allsp, k)))
    return Set(splits for (((splits,),kwd),ret) in helper(ss) if ret)

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
