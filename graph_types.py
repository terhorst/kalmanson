from sage.all import graphs, parallel, permutations, Set
from itertools import combinations
from memoized import memoized

@memoized
def graph_types(n, k):
    ret = []
    edges = filter(lambda (x,y): (x-y) % n not in [1, n-1], 
                combinations(range(n), 2))
    for e in combinations(edges, k):
        g = graphs.CycleGraph(n)
        g.add_edges((a,b,'split') for a,b in e)
        if all(not h.is_isomorphic(g) for h in ret):
            ret.append(g)
    return ret

def match_split_system(n, ss):
    perms = permutations(range(n))
    k = len(ss)
    gt = graph_types(n,k)
    return [g for g in gt if
            any(permute_ss(p, ss) == Set(edge_to_split(n, e) for e in split_edges(g))
                for p in perms)]

def permute_ss(p, ss):
    n = len(p)
    m = dict(zip(range(n), p))
    return Set(Set(Set(map(m.get, A)) for A in sp) for sp in ss)

def split_edges(g):
    return [(a,b) for a,b,c in g.edges() if c=="split"]

@memoized
def edge_to_split(n, edge):
    a,b = edge
    assert a < b
    A = Set(range(a,b))
    return Set([A, Set(range(n)) - A])
