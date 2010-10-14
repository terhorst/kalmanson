import sys, os, re
# import multiprocessing
# import multiprocessing.dummy as multiprocessing
from multiprocessing_jth import pyprocessing
from matrify import matrify
import numpy as np
from itertools import izip_longest, combinations, ifilter, imap
import itertools as it
from functools import partial
from sage.all import matrix, binomial, zero_matrix, combinations_iterator, vector
from memoized import memoized
import utility_matrices as um

np.set_printoptions(linewidth=100)

# Utility functions
def collapse_list(lst):
    "Turn a list of numbers into a dict counting how many times each number occurs."
    ret = {}
    for num in lst:
        ret[num] = ret.get(num,0) + 1
    return ret

@matrify
def conjugate(D,M):
    "Return matrix conjugation (D^T)*M*D."
    return np.dot(np.array(D).T, np.dot(M, D))
    

def textual_vector(n):
    return textual_symmetric_matrix(n)[triu_indices(n,1)]

@matrify
def symmetric_matrix(n, triu):
    "Make a symmetric nxn matrix from upper triangular entries triu."
    A = np.zeros([n,n], dtype=int)
    A[triu_indices(n,1)] = triu
    return A + A.T

def rowsums(M):
    return map(sum, M)

def number_of_nonzero(vec):
    return np.nonzero(vec)[0].shape[0]

def invert_dict(d):
    return dict([v,k] for k,v in d.items())

def prettify_cone_descs(f):
    "Prettify the cone descriptions of a fan."
    rays_h = {}
    for i,r in enumerate(f.rays()):
        rays_h[r] = "r%i" % i
    for c in flatten(f.cones()):
        c.rename(",".join(map(lambda r: rays_h[r], c.rays())))
        
def distance_matrix(n):
    "Make a stock distance matrix."
    M = np.zeros([n,n], dtype=int)
    M[triu_indices(n, 1)] = np.arange(binomial(n,2))
    M += M.T
    return matrix(M)

def textual_symmetric_matrix(n):
    """
    Return a symmetric matrix whose (i,j) entry is "d_ij".
    """
    return np.frompyfunc(lambda x,y: "d_%i%i" % (x+1,y+1), 2, 1).outer(range(n), range(n))

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def triu_indices(n, d=0):
    """
    Return the indices for the upper triangular part of nxn
    square matrix with diagonal offset d. (Backport of a
    Numpy 1.5 function.)
    """
    r = []
    c = []
    for i in range(n):
        for j in range(i+d, n):
            r.append(i)
            c.append(j)
    return list(map(np.array, (r,c)))

def upper_triangle(M):
    "Return the upper triangle of M as a vector."
    M = np.array(M)
    return vector(list(M[triu_indices(M.shape[0], 1)]))

# Real code
class KalmansonSystem:
    def __init__(self, n, ieqs=None):
        self._n = n
        if ieqs:
            assert ieqs.ncols() == binomial(n,2)
            self._ieqs = ieqs
        else:
            self._ieqs = kalmanson_matrix(n)

    def __eq__(self, other):
        try:
            return self.to_Set() == other.to_Set()
        except:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def to_Set(self):
        """
        Return a set consisting of the rows of self. Needed to get
        around immutability issues.
        """
        rows = self._ieqs.rows()
        [row.set_immutable() for row in rows]
        return Set(rows)

    def permute(self, g):
        """
        Returns a new KalmansonSystem obtained by permuting the leaves
        of self according to permutation p.
        """
        ind = permutation_vector(g)
        ieqs = matrix(np.vstack([np.array(row)[ind] for row in self._ieqs.rows()]))
        return KalmansonSystem(self._n, ieqs)

    def stabilizer(self):
        """
        Return the set of permutations in S_n under whose action
        the set of inequalities defining this system is unchanged.
        """
        S_n = SymmetricGroup(self._n)
        return filter(lambda g: self.permute(g)==self, S_n)

    def __repr__(self):
        strs = textual_symmetric_matrix(self._n)[triu_indices(self._n, 1)]
        vars = vector(var(",".join(strs)))
        descs = [str(vars * row) + " >= 0" for row in self._ieqs]
        return "A Kalmanson system of inequalities:\n%s" % "\n".join(descs)

def alt_kalmanson_matrix(n):
    r,c = triu_indices(n, 1)
    ncols = len(r)
    pos = np.arange(ncols)
    row_ind = lambda (i,j): pos[np.logical_and(r==i-1, c==j-1)]
    upright = lambda (i,j): (i,j) if i<j else (j,i)
    get_ind = lambda ind_lst: np.array(map(row_ind, map(upright, ind_lst)))
    rows = []
    kal_range = lambda m: range(m+1, n+1)
    for i in range(1, n+1):
        for j in kal_range(i):
            for k in kal_range(j):
                for l in kal_range(k):
                    for inds in [((i,j), (k,l), (i,k), (j,l)),
                            ((i,l), (j,k), (i,k), (j,l))]:
                        indices = get_ind(inds)
                        row = np.zeros(ncols, dtype=np.int)
                        row[indices] = [-1, -1, 1, 1]
                        rows.append(row)

    return matrix(np.vstack(rows))

def kalmanson_matrix(n, aug=False):
    r,c = triu_indices(n, 1)
    k = len(r)
    inds = np.arange(k)
    row_ind = lambda (i,j): inds[np.logical_and(r==i-1, c==j-1)]
    upright = lambda (i,j): (i,j) if i<j else (j,i)
    get_ind = lambda ind_lst: np.array(map(row_ind, map(upright, grouper(2, ind_lst))))

    rows = []
    for i in range(1,n-2):
        for j in range(i+2, n):
            ind_lst = (i, j+1, i+1, j, i, j, i+1, j+1)
            indices = get_ind(ind_lst)
            row = np.zeros(k, dtype=np.int)
            row[indices] = [-1, -1, 1, 1]           
            rows.append(row)

    sub = len(rows)

    for i in range(2, n-1):
        ind_lst = (i, 1, i+1, n, i, n, i+1, 1)
        indices = get_ind(ind_lst)
        row = np.zeros(k, dtype=np.int)
        row[indices] = [-1, -1, 1, 1]           
        rows.append(row)

    mat = matrix(np.vstack(rows))
    mat.subdivide(sub, None)
    if aug:
        return zero_matrix(len(rows),1).augment(mat)
    else:
        return mat

def _check_ind(n):
    M = []
    for i in range(1, n-2):
        for j in range(i+2, n):
            M.append(A(n,i,j))
    return(M)

def permute_matrix(g, M):
    mat = permutation_action(g, permutation_action(g, M).transpose()).transpose()
    mat.set_immutable()
    return mat

def orbit(M):
    "Return the orbit of the nxn symmetric matrix M."
    n = M.ncols()
    Sn = SymmetricGroup(n)
    d = {}
    for g in Sn:
        mat = permute_matrix(g, M)
        try:
            d[mat].append(g)
        except KeyError:
            d[mat] = [g]
    return d

def stabilizer(M):
    "Return the stabilizer subgroup of the nxn symmetric matrix M."
    n = M.ncols()
    Sn = SymmetricGroup(n)
    return filter(lambda g: permute_matrix(g, M)==M, Sn)

def permutation_vector(g):
    """
    The distance matrix obtained by permuting taxa according
    to the permutation g \in SymmetricGroup(n).
    """
    n = g.parent().degree()
    M = distance_matrix(n)
    M = permutation_action(g,permutation_action(g,M).transpose()).transpose()
    perm_vec = list(np.array(M)[triu_indices(n, 1)])
    return perm_vec
        
def permuted_kalmanson_matrix(A, perm):
    pv = permutation_vector(perm)
    M = matrix(np.vstack([np.array(r)[pv] for r in A.rows()]))
    M.set_immutable()
    return M
   
def set_of_kalmanson_matrices(n):
    "Return the set of all distinct Kalmanson nxn matrices."
    A = kalmanson_matrix(n)
    Sn = SymmetricGroup(n)
    mats = [permuted_kalmanson_matrix(A,p) for p in Sn]
    Kn = Set(map(set_of_rows, mats))
    return Kn

def kalmanson_polyhedra(n):
    """
    Return a list of all Kalmanson polyhedra for n taxa (including reorderings.)
    """
    R = rays(n)
    return [Polyhedron(rays=map(upper_triangle, map(partial(permute_matrix, g), R))) \
                for g in SymmetricGroup(n)]

def make_cone(p):
    from kalmanson import upper_triangle
    return Cone(map(upper_triangle, p), ZZ**binomial(p[0].ncols(), 2))

def kalmanson_cones(n):
    R = rays(n)
    ray_sets = Set([Set([permute_matrix(g,M) for M in R]) for g in SymmetricGroup(n)])
    p_iter = sage.parallel.use_fork.p_iter_fork(sage.parallel.ncpus.ncpus() * 2,30)
    P = parallel(p_iter=p_iter)
    return [ret for ((poly,kwd),ret) in P(make_cone)(ray_sets.list())]

def kalmanson_fan(n):
    return Fan(Set(kalmanson_cones(n)))

# @parallel(ncpus=2, p_iter="multiprocessing")
def std_kalmanson_polyhedron(n, lineality=False):
    return kalmanson_polyhedron(kalmanson_matrix(n), lineality)

def rays_lines_M2(M):
    macaulay2.eval('load "Polyhedra.m2"')
    P = macaulay2.matrix(M).intersection()
    rays,lines = [m.to_sage().transpose() for m in [P.rays(), P.linSpace()]]
    return rays,lines

def kalmanson_polyhedron(M, lineality=False):
    rays,lines = [M.rows() for M in rays_lines_M2(M)]
    if lineality:
        return Polyhedron(rays=rays, lines=lines)# , eqns=ortho)
    else:
        return Polyhedron(rays=rays)

def projected_distances(M, v):
    """
    Return the distance vector after modding out by the lineality space
    of Kalmanson polyhedron given by the Kalmanson inequality matrix M.
    """
    L = M.right_kernel().basis_matrix()
    B = M.stack(L).transpose()
    P = matrix(np.diag([1]*M.nrows() + [0]*L.nrows()))
    return B*P*B**(-1)*v

def positive_lineality_space(n):
    A = kalmanson_matrix(n)
    n = A.nrows()
    c = A.ncols()
    eq = zero_matrix(n,1).augment(A).rows()
    ieq = zero_matrix(c,1).augment(identity_matrix(c)).rows()
    return Polyhedron(eqns=eq, ieqs=ieq)

def rays(n):
    "The rays associated with the standard Kalmanson inequalities on n taxa."
    ret = []
    for i in range(2, n-1):
        ret.append(um.V(n, i))

    for i in range(1, n-2):
        for j in range(i+2, n):
            ret.append(um.V(n,i,j))

    return ret

def textual_rays(n):
    tsm = textual_symmetric_matrix(5)
    return map(list, [tsm[np.nonzero(np.triu(r,1))] for r in rays(n)])

def block_structure(M):
    "Diagonal block structure of nxn symmetric matrix M."
    lst = []
    n = M.ncols()
    i = 0
    while i < n:
        subs = [M.submatrix(i,i,j,j) for j in range(1, n-i+1)]
        block_size = [S for S in subs if S.is_zero()][-1].ncols()
        lst.append(block_size)
        i += block_size
    return tuple(lst)

