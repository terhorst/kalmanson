import sys, os, re
# import multiprocessing
# import multiprocessing.dummy as multiprocessing
from multiprocessing_jth import pyprocessing
from matrify import matrify
import numpy as np
from itertools import izip_longest, combinations, ifilter, imap, product
import itertools as it
from functools import partial
from sage.all import matrix, binomial, zero_matrix, \
       combinations_iterator, vector, permutation_action, db, Set, Fan
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

def orbit(n, M):
    "Return the orbit of the nxn symmetric matrix M."
    Sn = SymmetricGroup(n)
    return Set([permute_action(g, M) for g in Sn])

def standardize_elt(v):
    v = vector(v)
    return v if v[0]==1 else v*-1

def orbit_by_element(g,v):
    v = tuple(v)
    initial = copy(v)
    v = tuple(standardize_elt(g(v)))
    orb = []
    while initial != v:
        orb.append(copy(v))
        v = tuple(standardize_elt(g(v)))
    orb.append(copy(v))
    return Set(orb)

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
    Return a list of all Kalmanson polyhedra for n taxa (including
    reorderings.)
    """
    R = rays(n)
    return [Polyhedron(rays=map(upper_triangle, map(partial(permute_matrix, g), R))) \
                for g in SymmetricGroup(n)]

def make_cone(p):
    from kalmanson import upper_triangle
    return Cone(map(upper_triangle, p), ZZ**binomial(p[0].ncols(), 2))

def kalmanson_cones(n):
    desc = "kalmanson_cones_%i" % n
    try:
        return db(desc)
    except IOError:
        R = rays(n)
        ray_sets = Set([Set([permute_matrix(g,M) for M in R]) for g in SymmetricGroup(n)])
        p_iter = sage.parallel.use_fork.p_iter_fork(sage.parallel.ncpus.ncpus() * 2,30)
        P = parallel(p_iter=p_iter)
        cones = [ret for ((poly,kwd),ret) in P(make_cone)(ray_sets.list())]
        db_save(Set(cones), desc)
        return cones

def kalmanson_fan(n):
    fn = "kalmanson_fan_%i" % n
    try:
        return db(fn)
    except IOError:
        f = Fan(Set(kalmanson_cones(n)))
        f.db(fn)
        return f

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

def all_rays(n, Cn=None):
    if not Cn:
        Cn = kalmanson_cones(n)
    return Set([Set([tuple(ray_sign_vector(symmetric_matrix(n, r))) for r in cone.rays()]) \
            for cone in Cn])

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

def ray_sign_pattern(R):
    d = {-1:"-", 1:"+"}
    sv = ray_sign_vector(R)
    return "".join(map(d.get, sv))

def ray_sign_vector(R):
    def SignAlternator():
        i = -1
        while True:
            i *= -1
            yield i
    sa = SignAlternator()
    v = vector([s for s,i in zip(sa, block_structure(R)) for k in range(i) ])
    v.set_immutable()
    return v

def test_permutation(g,R):
    Rp = permute_matrix(g,R)
    pat = map(ray_sign_pattern, [R,Rp])
    vec = np.array(permutation_action(g, ray_sign_vector(R)))
    if vec[0]==-1:
        vec *= -1
    pred = vector(list(vec))
    d = {-1:"-", 1:"+"}
    pred = "".join(map(d.get, pred))
    #print R,"\n\n",Rp
    #print "p0: %s\tp': %s\tp_pred: %s" % (pat[0], pat[1], pred)
    assert pred==pat[1]

def cone_sign_matrix(n, C):
    rays = [symmetric_matrix(n, r) for r in C.rays()]
    return matrix([ray_sign_vector(r) for r in rays])

def ray_graph(R):
    G = Graph(R)
    n = R.ncols()
    # splits = np.split(np.arange(R.ncols(), dtype=int), np.cumsum(block_structure(R)))
    # partitions = [list(sp) for sp in splits if list(sp)]
    pluses = tuple(i for i,c in enumerate(ray_sign_vector(R)) if c==1)
    partitions = [pluses, tuple(i for i in range(n) if i not in pluses)]
    return G,partitions

def ray_split(R):
    sv = ray_sign_vector(R)
    pluses = tuple(i for i,c in enumerate(sv) if c==1)
    partition = [pluses, tuple(i for i in range(len(sv)) if i not in pluses)]
    return partition

def permutations_which_fix(n,m):
    Rn = Set(rays(n))
    return filter(lambda g: Set(permute_matrix(g,M) for M in Rn).intersection(Rn).cardinality() == m,
            SymmetricGroup(n))
            
def pairwise_compatible(sp_A, sp_B):
    "True if the two split sets sp_A and sp_B are pairwise compatible."
    return any(Set(A).intersection(Set(B)).cardinality() == 0 \
            for A,B in it.product(sp_A, sp_B))

def compatibility_degree(splits):
    "Check that this list of splits is pairwise compatible."
    return sum(1 for pair in it.combinations(splits,2) \
            if not any(Set(A).intersection(Set(B)).cardinality() == 0
                for A,B in it.product(*pair)))

def cone_splits(C,n):
    return [G[1] for G in make_cone_graph(C,n)]

def make_cone_graph(C,n):
    graphpairs = [ray_graph(symmetric_matrix(n,r)) for r in C.rays()]
    # return [G.plot(partition=part, layout="circular") for G,part in graphpairs]
    return graphpairs

def ray_magic(C,n):
    ray_sets = [Set([ray_sign_vector(symmetric_matrix(n,r)) for r in cone]) for cone in C] 
    canonical = ray_sets[0]
    d = {}
    for g in SymmetricGroup(n):
        permute = Set([tuple(permutation_action(g, v)) for v in canonical])
        d[permute] = d.get(permute, []) + [g]
    return d

def ray_matrices(cone,n):
    return [symmetric_matrix(n, r) for r in cone.rays()]

def all_ray_matrices(n):
    Cn = kalmanson_cones(n)
    return uniq(M for cone in Cn for M in ray_matrices(cone, n))

def graph_rays_cones(C):
    n = C[0].lattice_dim() / 2
    Rn = Set([r for g in SymmetricGroup(n) for v in map(ray_sign_vector, rays(n))
            for r in orbit_by_element(g,v)])
    d = {-1: "-", 1: "+"}
    signer = lambda v: "".join(map(d.get, v))
    descs = Set(map(signer, Rn))
    part = [descs, []]
    G = Graph()
    G.add_vertices(descs)
    i = 0
    for c in C:
        i += 1
        desc = "Cone %i" % i 
        part[1].append(desc)
        G.add_vertex(desc)
        for r in c.rays():
            sv = ray_sign_pattern(symmetric_matrix(n, r)) 
            G.add_edge(sv, desc)
    
    return G,part

def ray_to_splits(n, ray):
    mat = symmetric_matrix(n, ray)
    vec = ray_sign_vector(mat)
    ns = range(n)
    pos = Set(ind for s,ind in zip(vec, ns) if s==1)
    neg = Set(ns) - pos
    A = pos if 0 in pos else neg
    return Set([A, Set(ns) - A])
    if pos.cardinality() < neg.cardinality():
         return pos
    elif pos.cardinality() == neg.cardinality():
         return pos if sorted(pos)[0] < sorted(neg)[0] else neg
    else:
         return neg

def non_trivial_splits(n):
    rng = Set(range(1,n+1))
    splits = Set([Set(tuple(s)) for s in powerset(rng) if len(s)>1 and len(s)<n-1])
    splits = Set([x if x.cardinality() <= floor(n/2) else rng - x \
        for x in splits])
    splits = Set([rng - x if x.cardinality() == n/2 and min(x) > min(rng - x) else x
        for x in splits])
    return splits


def show_partition_types(n,k):
    G = []
    Cn = graphs.CycleGraph(n)
    for edges in combinations_iterator(filter(lambda (x,y): abs(x-y) % (n-1) > 1, combinations_iterator(range(n), 2)), k):
        g = Cn.copy()
        g.add_edges(edges)
        G.append(g)
    H = []
    for g in G:
        match = False
        for j in H:
            if g.is_isomorphic(j):
                match = True
                break
        if not match:
            H.append(g)
    return H

def f_vector(fan):
    return map(len, fan.cone_lattice().level_sets())[:-1]

def h_vector(fan):
    f = f_vector(fan)
    d = len(f)
    x = var('x')
    poly = sum([f[i]*(x-1)**(d-i) for i in range(d)])
    return [a for a,b in reversed(poly.expand().coefficients())]

def number_simplex(n, rays, sought=None):
    # General recurrence relation for multinomial:
    # multi(k1,...,kn) = multi(k1-1,k2,...,kn) + multi(k1,k2-1,...,kn) + ... + multi(k1,...,kn-1)

    dim = len(rays)

    def level_set(c):
        if len(c) == dim:
            n = _simplicial_number(c,rays)
            if n==sought:
                print Exception("Found %i: f(%s,%s)" % (sought, c, rays))
            return n
        else:
            subsimplices = [(c[:-1] + (c[-1]-i,) + (i,)) for i in range(c[-1]+1)]
            return map(level_set, subsimplices)

    return [level_set((i,)) for i in range(n)]

@memoized
def _simplicial_number(coords,rays):
    z = list(np.nonzero(coords)[0])
    if len(z) == 0:
        return 1
    elif len(z) == 1:
        return rays[z[0]]
    elif len(z) < len(rays):
        z = coords.index(0)
        cp = coords[:z] + coords[z+1:]
        rp = rays[:z] + rays[z+1:]
        return f(cp, rp)
    else:
        dim = len(coords)
        shift_vec = np.array([-1] + [0]*(dim-1))
        return sum(f(tuple(np.roll(shift_vec, i) + coords),rays) for i in range(dim))

def enumerate_splits(n,k,invalid=False):
    f = invalid_splits if invalid_splits else valid_splits
    s = f(n,k)
    tups = [tuple(sorted(map(len, x))) for x in s]
    return dict((u,tups.count(u)) for u in uniq(tups))

def valid_splits(n,k):
    faces = kalmanson_fan(n)(k)
    sp = Set([Set([ray_to_splits(n, r) for r in f.rays()]) for f in faces])
    return sp

def invalid_splits(n,k):
    faces = kalmanson_fan(n)(k)
    return Set(map(Set, combinations(non_trivial_splits(n),k))) - \
            valid_splits(n,k)

def dimension_of_span(n, k):
    dim_counts = {}
    for rays in combinations_iterator(all_ray_matrices(n), k):
        kers = [r.image() for r in rays]
        dim = reduce(lambda x,y: x.intersection(y), kers).dimension()
        dim_counts[dim] = dim_counts.get(dim,0) + 1
    return dim_counts

def count_stuff(fan):
    cones = fan.cones()
    d = len(cones)
    for i in range(1,d):
        c_i = cones[i]
        print "Dimension %i (%i cones):" % (i, len(c_i))
        dims = {}
        for c in c_i:
            rays = c.rays()
            ind = tuple(sum(1 for cone in cones[j] if all(ray in cone.rays() for ray in rays)) \
                    for j in range(i, d))
            dims[ind] = dims.get(ind,0) + 1
        print "\n".join("\t%s: %i" % (k,v) for k,v in dims.iteritems())

def cones_fun(n, k):
    fan = kalmanson_fan(n)
    rays = [c.rays()[0] for c in fan(1)]
    cones = Set(fan(k))
    for r in combinations(rays, k):
        c = Cone(r) 
        if c in cones:
            print prod(symmetric_matrix(n, ray) for ray in r).change_ring(GF(2))
            print ""

def cone_to_split_system(n, cone):
    return Set(ray_to_splits(n, r) for r in cone.rays())
