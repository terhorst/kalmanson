from sage.all import graphs, parallel, permutations, Set, floor
from itertools import combinations
from memoized import memoized
from kalmanson import kalmanson_fan, cone_to_split_system

@memoized
def graph_types(n, k):
    ret = []
    # edges = combinations(range(n), 2)
    edges = [[i, j] for i,j in combinations(range(n), 2)
                if abs(i-j) != 1 and abs(i-j) != (n-1)]
    for e in combinations(edges, k):
        g = graphs.CycleGraph(n)
        g.add_edges((a,b,'split') for a,b in e)
        if all([not h.is_isomorphic(g) for h in ret]):
            ret.append(g)
    return ret

@parallel(p_iter='multiprocessing')
def match_split_system(ss):
    from graph_types import graph_types, permute_ss, split_edges, edge_to_split
    n = sum(map(len, ss[0]))
    perms = permutations(range(n))
    k = len(ss)
    gt = graph_types(n,k)
    return [g for g in gt if
            any([permute_ss(p, ss) ==
                Set([edge_to_split(n, e)
                    for e in split_edges(g)])
                for p in perms])]

def permute_ss(p, ss):
    n = len(p)
    m = dict(zip(range(n), p))
    return Set([
            Set(
                [Set(map(m.get, A)) for A in sp]
                )
            for sp in ss])

def split_edges(g):
    return [(a,b) for a,b,c in g.edges() if c=="split"]

@memoized
def edge_to_split(n, edge):
    a,b = edge
    assert a < b
    A = Set(range(a,b))
    return Set([A, Set(range(n)) - A])

def kalmanson_graph_types(n, k):
    cones = kalmanson_fan(n)(k)
    ss = [cone_to_split_system(n, cone) for cone in cones]
    return list(match_split_system(ss))
