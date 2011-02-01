from multiprocessing import Pool, cpu_count
from consecutive_ones import circular_ones
from matrify import matrify
import numpy as np
from memoized import memoized
from itertools import combinations, imap, izip, repeat, starmap, product
from sage.all import matrix, SymmetricGroup, permutation_action, uniq, \
        Set, partitions, multinomial, vector
from progressbar import ProgressBar

M1_COL = [(1,0,1),(0,1,1),(1,1,0)]
ALL_COL = list(product([0,1], repeat=3))
NONM1_COL = [c for c in ALL_COL if c not in M1_COL]

def __base_matrix(n,k):
    base = [1,1] + [0]*n
    return np.vstack([np.roll(base,i) for i in range(k)])

@matrify
def m1(n):
    return __base_matrix(n,n+2)

@matrify
def m2(n):
    base2 = [0] + [1]*(n+2)
    mat = np.hstack([__base_matrix(n,n+1), np.zeros([n+1,1],dtype=int)])
    return np.vstack([mat, [np.roll(base2,i) for i in [n+1,0]]])

@matrify
def m3(n):
    mat = np.hstack([__base_matrix(n,n+1), np.zeros([n+1,1],dtype=int)])
    return np.vstack([mat, ([0] + [1]*n + [0,1])])

@matrify
def m4():
    base = [1,1] + [0]*4
    return np.vstack([np.roll(base,i*2) for i in range(3)] + [[0,1]*3])

@matrify
def m5():
    return [[1,1,0,0,0], [1,1,1,1,0], [0,0,1,1,0], [1,0,0,1,1]]

@memoized
def split_to_vector(n, split):
    v = np.zeros(n, dtype=int)
    split = split if 0 not in split else \
        [i for i in range(n) if i not in split]
    v[list(split)]=1
    return v

def all_submatrices(M, m, n):
    return (M.matrix_from_rows_and_columns(list(r),list(c))
            for r in combinations(range(M.nrows()), m)
            for c in combinations(range(M.ncols()), n))

def matrix_contains(M, S):
    "Does M contain a configuration of S as a submatrix?"
    m,n = (S.nrows(), S.ncols())
    configs = submatrix_configurations(S)
    return any(sm in configs for sm in all_submatrices(M, m, n))

def orbit(S):
    m,n = (S.nrows(), S.ncols())
    Sr = SymmetricGroup(m)
    Sc = SymmetricGroup(n)
    return [ 
        permutation_action(h, permutation_action(g, S).transpose()).transpose()
        for g in Sr for h in Sc
        ]

def col_permutation_action(g, M):
    M = permutation_action(g, M.transpose()).transpose()
    M.set_immutable()
    return M

def row_orbit(M):
    mats = [permutation_action(g, M) for g in SymmetricGroup(M.nrows())]
    [m.set_immutable() for m in mats]
    return mats

@memoized
def submatrix_configurations(S):
    configs = orbit(S)
    [m.set_immutable() for m in configs]
    return uniq(configs)

@memoized
def ss_to_matrix(n, ss):
    M = matrix(np.vstack([split_to_vector(n,sp) for sp in ss]))
    M.set_immutable()
    return M

def __incomp_helper(tup):
    n,M,S = tup
    return (M,matrix_contains(M,S))

def incomp_tucker(n,k,S,parallel=True):
    # ss = ProgressBar()(circular_ones(n,k,True))
    ss = [ ss_to_matrix(n,ss) for ss in circular_ones(n,k,True) ]
    if parallel:
        p = Pool(cpu_count())
        args = izip(repeat(n),ss,repeat(S))
        iter = p.imap_unordered(__incomp_helper, args, 1000)
    else:
        args = izip(repeat(n),ProgressBar()(ss),repeat(S))
        iter = imap(__incomp_helper, args)
    return [ss for ss,cont in iter if cont]

def matrices_with_cols(mats, cols):
    cols = map(tuple, cols)
    return filter(lambda M: sorted(cols) in
            imap(sorted, combinations(map(tuple, M.columns()), len(cols))),
            mats)

def interesting_permutations(S):
    m,n = (S.nrows(), S.ncols())
    Sr = SymmetricGroup(m)
    Sc = SymmetricGroup(n)
    try:
        ( (g,h) for g in Sr for h in Sc 
            if g != Sr.identity() and h != Sc.identity() and 
            S == permutation_action(h, permutation_action(g, S).transpose()).transpose() ).next()
    except StopIteration:
        return 0
    return 1

@memoized
def classify_mats(mats):
    d = {}
    numones = {}
    sys_col = [(1,0,1),(0,1,1),(1,1,0)]
    pb = ProgressBar()(mats)
    for m in pb:
        cols = map(tuple, list(m.columns())[:-1])
        map(cols.remove, sys_col)
        M = matrix(cols).transpose()
        M.set_immutable()
        Mcs = submatrix_configurations(M)
        try:
            k = (k for k in d.keys() if k in Mcs).next()
            d[k] += 1
        except StopIteration:
            d[M] = 1
    return d

def mat_class_pp(d):
    n = d.keys()[0].ncols() + 3
    mtable = dict(zip(starmap(multinomial, partitions(n)), partitions(n)))
    for k,v in sorted(d.iteritems(), key=lambda tup:tup[1]):
        cols = k.columns()
        rows = k.rows()
        [c.set_immutable() for c in cols]
        [r.set_immutable() for r in rows]
        print "%s %s: %i \n" % (k, mtable.get(v, "??"), v)

def class_mats_m1(lst):
    d = {0:[],1:[],2:[]}

    for M in lst:
        x = sum(min(1, M.columns().count(vector(c))) for c in SYS_COL)
        d[x].append(M)

    return d

def construct_lift(M, v):
    M2 = M.matrix_from_columns(range(M.ncols() - 1))
    M2 = M2.augment(matrix([v, [0]*3]).transpose())
    M2.set_immutable()
    return M2

