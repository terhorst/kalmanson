from sage.all import Set, combinations_iterator, parallel, binomial
import itertools as it

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

    def merge(self, other_split):
        return Split(self.X, self.A.intersection(other_split.A)) 

    def blocks(self):
        return Set([self.A, self.X - self.A])

class SplitSystem(object):
    def __init__(self,splits=[]):
        self._splits = Set(splits)

    def add_split(self, split):
        self._splits += split

    def splits(self):
        return self._splits

    def is_compatible(self):
        return all(
                any(A.intersection(B).cardinality()==0 \
                for A,B in self._block_product(double)
                for double in self._split_combinations(2)))

    def is_weakly_compatible(self):
        return all(
                A1.intersection(A2).intersection(A3).cardinality() == 0 or
                (A1 - A2 - A3).cardinality() == 0 or
                (A2 - A1 - A3).cardinality() == 0 or
                (A3 - A1 - A2).cardinality() == 0
                for triple in self._split_combinations(3)
                for A1,A2,A3 in self._block_product(triple))

    def _block_product(self, lst):
        return it.product(*[x.blocks() for x in lst])

    def _split_combinations(self, n):
        return combinations_iterator(self._splits, n)

    def is_circular(self):
        return SplitSystem(self._splits.union(Set(S1.merge(S2)
                    for S1,S2 in self._split_combinations(2)))).is_weakly_compatible()

def all_splits(X):
    X = Set(X)
    return Set(Split(X,b) for i in range(2, X.cardinality() - 1) \
            for b in combinations_iterator(X, i))

@parallel('multiprocessing')
def __circular_splits_helper(sp): 
    from splits import SplitSystem
    return SplitSystem(sp).is_circular()

@parallel('multiprocessing')
def __weak_splits_helper(sp):
    from splits import SplitSystem
    return SplitSystem(sp).is_weakly_compatible()

def weakly_compatible_splits(n, k):
    return __splits_checker(n, k, __weak_splits_helper)

def circular_splits(n, k):
    return __splits_checker(n, k, __circular_splits_helper)

def __splits_checker(n, k, helper):
    S = all_splits(range(1, n+1))
    J = binomial(S.cardinality(),k)
    return [splits for (((splits,),kwd),ret)
            in helper(list(combinations_iterator(S, k))) if ret]
