from sage.all import graphs, parallel, permutations, Set, floor, SymmetricGroup, DihedralGroup
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
        if not any(is_equivalent(g,h) for h in ret):
            ret.append(g)
    return ret

@parallel(p_iter='multiprocessing')
def match_split_system(ss):
    from graph_types import graph_types, permute_ss, split_edges, edge_to_split
    n = sum(map(len, ss[0]))
    k = len(ss)
    gt = graph_types(n,k)
    perms = SymmetricGroup(n)
    ss_sig = sorted(map(split_len, ss))
    for g in gt:
        gss = graph_to_ss(g)
        gss_sig = sorted(map(split_len, gss))
        if gss_sig != ss_sig:
            continue
        elif any([permute_ss(p, ss) == gss for p in perms]):
            return g
    # return [g for g in gt if any([permute_ss(p, ss) == Set([edge_to_split(n, e)
    #                for e in split_edges(g)]) for p in perms])]

def split_len(sp):
    return min(map(len, sp))

def ss_signature(ss):
    return sorted(map(split_len, ss))

def is_equivalent(g,h):
    gss,hss = map(graph_to_ss, [g,h])
    n = len(g.vertices())
    return any(permute_ss(p,gss) == hss for p in DihedralGroup(n))

@memoized
def graph_to_ss(g):
    n = len(g.vertices())
    return Set([edge_to_split(n,e) for e in split_edges(g)])

def permute_ss(p, ss):
    return Set([
            Set(
                [Set(map(lambda x: p(x+1) - 1, A)) for A in sp]
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
