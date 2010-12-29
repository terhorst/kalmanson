from sage.all import graphs
from itertools import combinations

def all_graph_types(n, k):
    edges = filter(lambda (x,y): (x-y) % n not in [1, n-1], combinations(range(n), 2))
    ret = []
    for e in combinations(edges, k):
        g = graphs.CycleGraph(n)
        g.add_edges(e)
        if all(not g.is_isomorphic(h) for h in ret):
            ret.append(g)
    return ret


