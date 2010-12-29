from sage.all import graphs, parallel
from itertools import combinations

def graph_types(n, k):
    ret = []
    edges = filter(lambda (x,y): (x-y) % n not in [1, n-1], combinations(range(n), 2))
    for e in combinations(edges, k):
        g = graphs.CycleGraph(n)
        g.add_edges(e)
        if all(not h.is_isomorphic(g) for h in ret):
            ret.append(g)
    return ret
