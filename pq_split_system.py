from pqtree import PQTree
from sage.all import Set
import split_system

class PQSplitSystem(split_system.SplitSystem):
    def is_circular(self):
        X = self._splits[0].X
        Smax = Set([max(X)+1])
        pqt = PQTree(X.union(Smax))
        sps = [sp.A for sp in self._splits]
        return pqt.ReduceAll(sps)
