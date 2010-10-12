import sys, os, re
# import multiprocessing
# import multiprocessing.dummy as multiprocessing
from multiprocessing_jth import pyprocessing
from matrify import matrify
import numpy as np
from itertools import izip_longest, combinations, ifilter, imap
import itertools as it
from functools import partial
from sage.all import matrix, binomial, zero_matrix, combinations_iterator
from memoized import memoized
import utility_matrices as um

np.set_printoptions(linewidth=100)

# Utility functions
@matrify
def conjugate(D,M):
    "Return matrix conjugation (D^T)*M*D."
    return np.dot(np.array(D).T, np.dot(M, D))
    

def textual_vector(n):
    return textual_symmetric_matrix(n)[triu_indices(n,1)]

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

def make_basis(n, flat=True):
    triu = triu_indices(n)
    E_b = []
    E_b += [E(n,i+1) for i in range(n)]
    E_b += [E(n,i+1,i+1) for i in range(n)]
    if flat:
        return matrix(np.vstack([A.flatten() for A in E_b]))
    else:
        return E_b

def make_Aflat(n):
    M = [A(n,i) for i in range(2, n-1)] 
    for i in range(1, n-2):
        for j in range(i+2, n):
            M.append(A(n,i,j))
    triu = triu_indices(n,1)
    M = [T[triu].flatten() for T in M]
    return matrix(np.vstack(M))

def reformat_matrix(M, leading="0"):
    return re.sub('\\[\s*([^\\]]*)\\]', '%s \\1' % leading, M.str())

def polymake_kalmanson(n, ineq=False, fn=None):
    txt = ""
    if ineq:
        M = kalmanson_matrix(n,1)
        ortho = matrix(M).right_kernel().basis_matrix()
#        M = matrix(np.vstack([M, ortho, -ortho]))
        txt = "INEQUALITIES\n%s\nEQUATIONS\n%s"
        args = map(reformat_matrix, [matrix(M), ortho])
        out = txt % tuple(args)
    else:
        M = make_Aflat(n)
        txt = "POINTS\n%s"

    if fn:
        with open(fn, 'w+') as f:
            f.write(out)
    else:
        print out

def cdd_kalmanson(n):
    M = kalmanson_system(n,1)
    ortho = matrix(M).right_kernel().basis_matrix()
    M = np.vstack([-M, ortho, -ortho])
    txt = """H-representation
begin
%i %i integer
%s
end"""
    print txt % (M.shape[0], M.shape[1] + 1, reformat_matrix(matrix(M)))

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
    # return [kalmanson_polyhedron(permuted_kalmanson_matrix(n, p)) \
    # for p in SymmetricGroup(n)]
    M = kalmanson_matrix(n)
    mats = []
    for perm in SymmetricGroup(n):
        perm_vec = permutation_vector(perm)
        mats.append(matrix(np.vstack([np.array(r)[perm_vec] for r in M.rows()])))
    return map(kalmanson_polyhedron, mats)
    # return [kalmanson_polyhedron(m) for p in kalmanson_polyhedron(mats)]
    # return kalmanson_polyhedron(mats)

def kalmanson_cones(n):
    polys = kalmanson_polyhedra(n)
    cones = [Cone(p, ZZ**(binomial(n, 2))) for p in polys]
    return cones

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

def ray_matrices(n):
    A,B = positive_lineality_space(n),std_kalmanson_polyhedron(n)
    return map(matrix, [m.rays() for m in (A,B)]) 

def build_ray_matrix(M):
    # There will by (nC2 - n) rays. Each is nC2-dimensional.
    n = M.nrows()
    m = M.ncols()
    A = zero_matrix(n, m)
    # The matrix contains n columns with only a -1.
    
    single_cols = lonely_columns(M)
    assert len(single_cols) >= A.nrows()
    single_cols = single_cols[:A.nrows()]
    
    cols = M.columns()

    for col in single_cols:
        row = np.nonzero(cols[col])[0][0]
        assert M[row,col] == -1
        A[row,col] = 1
        
    return A

# @parallel(ncpus=8, p_iter="multiprocessing")
def check_valid(A,Mnp,inds):
    inds = np.array(inds)[:,:,0]
    B = np.array(A, dtype=np.int)
    B[inds[:,0], inds[:,1]] = 1
    print B
    C = np.dot(Mnp,B.T)
    cond = np.all(C>=0) and np.all(np.sum(C,0) == 1) and np.all(np.sum(C,1)==1)
    return cond

def nonzero_entries(v):
    "Return the number of non-zero entries in the vector v."
    return nonzero_indices(v).shape[0]

def nonzero_indices(v):
    "Return the indices where v is nonzero."
    v = np.array(v)
    return np.nonzero(v)[0]

def lonely_columns(M):
    "Return column indices of M which have only one non-zero entry."
    return filter(lambda i: nonzero_entries(M[:,i])==1, range(M.ncols()))

def social_columns(M):
    return list(Set(range(M.ncols())) - Set(lonely_columns(M)))

def make_ones_matrix(dim, spots):
    r = len(spots)
    M = np.zeros([r, dim], dtype=int)
    for i,inds in enumerate(spots):
        M[i,inds] = 1
    return matrix(M).transpose()

def submatrices_with_column(A, col, k):
    rng = range(A.ncols())
    rng.remove(col)
    submat = A.matrix_from_columns(rng)
    return filter(partial(check_cols, A), \
            imap(lambda (x,y): Set((x,)+y), product([col], combinations(rng,k))))

def check_cols(A,cols):
    from kalmanson import rowsums, number_of_nonzero
    cols = list(cols)
    rs = rowsums(A.matrix_from_columns(cols))
    return sum(rs)==1 and number_of_nonzero(rs)==1

def pool_filter(pool, iter):
    a,b = it.tee(iter)
    return [cols for (mat,cols), keep in it.izip(a, pool.imap(check_cols, b)) if keep]

def look_for_patterns(n,m):
    A = kalmanson_matrix(n)
    tv = textual_vector(n)
    inds = range(A.ncols())

    p_iter = pyprocessing(8)
    P = parallel(p_iter=p_iter)
    args = list(it.product([A], combinations_iterator(range(A.ncols()), m)))
    cols = [vec for (((mat,vec),kwd),ret) in P(check_cols)(args) if ret]
    # cols = pool_filter(p, it.product([A], combinations_iterator(range(A.ncols()), m))) 
    # cols = filter(check_cols, combinations(range(A.ncols()), m))
    
    ret = []
    for col in cols:
        desc = tv[col]
        vec = np.zeros(binomial(n,2), dtype=int)
        vec[col] = 1
        ret.append([col, desc, vec])
        print "%s\t%s\t%s" % (desc, np.dot(A,vec), vec)
        print symmetric_matrix(n, vec)
    return ret 

def pretty_print_patterns(n):
    vecs = look_for_patterns(n)
    for k in sorted(vecs.keys()):
        print "%i\t%s\t%s\t%s\t%s" % (k, vecs[k][0][0], vecs[k][0][1], \
                vecs[k][1][0], vecs[k][1][1])
    return vecs

def ray_basis_from_paper(n):
    return {'i': [um.V(n, i+1) for i in range(2,n-1)],
            'ij': reduce(operator.add, \
                [[um.V(n, i, j) for j in range(i+2, n)] for i in range(1, n-2)])}

def textual_ray_basis(n):
    rb = ray_basis_from_paper(n)
    tm = textual_symmetric_matrix(n)
    vecs = dict([(k,[tm[np.nonzero(b)] for b in basis]) for (k,basis) in rb.iteritems()])
    return vecs

def ray_basis_pp(vecs, do_sort=True):
    if do_sort:
        print "\n".join(sorted(map(str, vecs)))
            #lambda a,b: cmp(a,b) if len(a)==len(b) else cmp(len(a),len(b))))
    else:
        print "\n".join(map(str,vecs))

@matrify
def type1_ray(n,m):
    M = np.array(um.E(n,m))
    M[m,:] = M[:,m] = 0
    return M

@matrify
def ray_D(n, m):
    return np.diag([1]*m + [0] + [1]*(n-m-1))

