from sage.all import Set, combinations_iterator, parallel, binomial
import itertools as it
from memoized import memoized

class Split(object):
    def __init__(self, X, A):
        A = Set(A)
        X = Set(X)
        assert A <= X
        self.X = X
        self.A = A if 1 in A else X - A
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

class SplitSystem(object):
    def __init__(self,splits=[]):
        self._splits = Set(splits)

    def __repr__(self):
        return "Split system: " + self._splits.__repr__()

    def __cmp__(self, other):
        return self._splits.__cmp__(other._splits)

    def __hash__(self):
        return self._splits.__hash__()

    def join(self, other):
        "Join self to the split system other."
        return SplitSystem(self._splits.union(other._splits))

    def add_split(self, split):
        self._splits += split

    def __len__(self):
        return self._splits.cardinality()

    def splits(self):
        return self._splits

    def is_compatible(self):
        return all(
                any(A.intersection(B).cardinality()==0 \
                for A,B in self._block_product(double)
                for double in self._split_combinations(2)))

    def is_weakly_compatible(self):
        return all(weakly_compatible_helper(Set(blocks))
                for triple in self._split_combinations(3)
                for blocks in self._block_product(triple))

    def _block_product(self, lst):
        return it.product(*[x.blocks() for x in lst])

    def _split_combinations(self, n):
        return combinations_iterator(self._splits, n)

    def is_circular(self):
        return SplitSystem(self._splits.union(Set(S1.merge(S2)
                    for S1,S2 in self._split_combinations(2)
                    if S1.B.intersection(S2.B).cardinality() > 0
                    ))).is_weakly_compatible()

def all_splits(n):
    X = Set(range(1, n+1))
    return Set(Split(X,b) for i in range(2, X.cardinality() - 1) \
            for b in combinations_iterator(X, i))

@memoized
def weakly_compatible_helper(S):
    A1,A2,A3 = S
    return A1.intersection(A2).intersection(A3).cardinality() == 0 or \
        (A1 - A2 - A3).cardinality() == 0 or \
        (A2 - A1 - A3).cardinality() == 0 or \
        (A3 - A1 - A2).cardinality() == 0

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

@memoized
def __splits_checker(n, k, helper):
    S = [SplitSystem([s]) for s in all_splits(n)]
    if k==1:
        return S
    Snm1 = __splits_checker(n, k-1, helper)
    ss = filter(lambda ss: len(ss)==k,
            (ss1.join(ss2) for ss1,ss2 in it.product(S, Snm1)))
    return Set(splits for (((splits,),kwd),ret) in helper(ss) if ret)
