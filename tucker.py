from consecutive_ones import circular_ones
from matrify import matrify
import numpy as np
from memoized import memoized
from itertools import combinations, imap
from sage.all import matrix, SymmetricGroup, permutation_action, uniq, Set

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

def submatrix_configurations(S):
    m,n = (S.nrows(), S.ncols())
    Sr = SymmetricGroup(m)
    Sc = SymmetricGroup(n)
    configs = [ permutation_action(h, permutation_action(g, S).transpose()).transpose()
            for g in Sr for h in Sc ]
    [m.set_immutable() for m in configs]
    return uniq(configs)

@memoized
def ss_to_matrix(n, ss):
    M = matrix(np.vstack([split_to_vector(n,sp) for sp in ss]))
    M.set_immutable()
    return M

def incomp_tucker(n,k,S):
    return filter(lambda M: matrix_contains(M,S),
            (ss_to_matrix(n,ss) for ss in circular_ones(n,k,True)))

def matrices_with_cols(mats, cols):
    cols = map(tuple, cols)
    return filter(lambda M: sorted(cols) in
            imap(sorted, combinations(map(tuple, M.columns()), len(cols))),
            mats)

def classify_mats(mats):
    d = {}
    for m in mats:
        cols = map(tuple, list(m.columns())[:-1])
        map(cols.remove, [(1,0,1),(0,1,1),(1,1,0)])
        M = matrix(sorted(matrix(cols).transpose()))
        M.set_immutable()
        d[M] = d.get(M,0) + 1

    for k,v in sorted(d.iteritems(), key=lambda tup:tup[1]):
        print "%s: %i \n" % (k,v)

    return d

