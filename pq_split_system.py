from pqtree import PQTree
from sage.all import Set
import split_system

class PQSplitSystem(split_system.SplitSystem):
    def is_circular(self):
        X = self._splits[0].X
        pqt = PQTree(X)
        sps = [sp.B for sp in self._splits]
        return pqt.ReduceAll(sps)
