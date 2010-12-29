from sage.all import graphs
from itertools import combinations

def all_graph_types(n, k):
    edges = filter(lambda (x,y): (x-y) % n not in [1, n-1], combinations(range(n), 2))
    ret = []
    for e in combinations(edges, k):
        g = graphs.CycleGraph(n)
        g.add_edges(e)
        match = False
        for h in ret:
            if h.is_isomorphic(g):
                match = True
                break
        if not match:
            ret.append(g)

    return ret


