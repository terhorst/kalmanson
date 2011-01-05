from sage.all import Set
from itertools import combinations, product
from memoized import memoized

@memoized
def weakly_compatible_helper(S):
    A1,A2,A3 = S
    return A1.intersection(A2).intersection(A3).cardinality() == 0 or \
        (A1 - A2 - A3).cardinality() == 0 or \
        (A2 - A1 - A3).cardinality() == 0 or \
        (A3 - A1 - A2).cardinality() == 0

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
        return SplitSystem(self._splits + Set([split]))

    def contains(self, split):
        return split in self._splits

    def __len__(self):
        return self._splits.cardinality()

    def __iter__(self):
        self._splits.__iter__()

    def splits(self):
        return self._splits

    def is_compatible(self):
        return all(
                any(A.intersection(B).cardinality()==0 \
                for A,B in self._block_product(double)
                for double in self._split_combinations(2)))

    def is_weakly_compatible(self):
        blockset = [Set(blocks)
            for triple in self._split_combinations(3)
            for blocks in self._block_product(triple)]
        return all(map(weakly_compatible_helper, blockset))

    def _block_product(self, lst):
        return product(*[x.blocks() for x in lst])

    def _split_combinations(self, n):
        return combinations(self._splits, n)

    def is_circular(self):
        ssp = SplitSystem(self._splits.union(Set(S1.merge(S2)
                    for S1,S2 in self._split_combinations(2)
                    if S1.B.intersection(S2.B).cardinality() > 0)))
        return ssp.is_weakly_compatible()
