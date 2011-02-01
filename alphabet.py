from memoized import memoized
from consecutive_ones import circular_ones
from tucker import ss_to_matrix, row_orbit
from sage.all import matrix, vector, prod, uniq
from itertools import permutations, product, imap
from progressbar import ProgressBar

ALPHA = dict(zip(((1,1,0),(1,0,1),(0,1,1),(1,0,0),
                 (0,1,0),(0,0,1),(0,0,0),(1,1,1)),
                 map(chr, range(ord("a"), ord("h") + 1))))

I = map(vector, [(1,1,0),(0,1,1),(1,0,1)])
II = map(vector, [(1,0,0),(0,1,0),(0,0,1)])
III = [vector((1,1,1))]

S = {
        'S1': (1, 0, 'X'),
        # 'S2': (1, 1, 0),
        'S3': (1, 0, 0),
        'S4': (1, 'X', 'X'),
        'S5': (0, 0, 'X'),
        # 'S6': (1, 1, 1),
        'S7': (0, 0, 0),
        'S8': (0, 'X', 'X'),
        # 'S9': (1, 1, 'X')
        }

def mat_to_word(ss):
    "".join(map(ALPHA.get, ss.columns()))

def __get_mats(n):
    return [ss_to_matrix(n, ss) for ss in circular_ones(n, 3)]

@memoized
def F(i,j,k,mats):
    return filter(lambda M: sum(vector(S) in M.columns() for S in I)==i and
                    sum(vector(S) in M.columns() for S in II)==j and
                    sum(vector(S) in M.columns() for S in III)==k,
                    mats)

def all_F(mats):
    d = dict(((i,j,k),F(i,j,k,mats)) for i in range(4) for j in range(5) for k in range(2))
    print_d(d)
    return d

def print_d(d):
    for k,v in d.iteritems():
        if len(v)>0:
            print "%s: %i" % (str(k),len(v))

def count_structure(n, sty="S1"):
    pat = S[sty]
    rows = []
    for i in pat:
        if i==1:
            rows.append([(0,) + ((1,)*(n-2)) + (0,)])
        elif i==0:
            rows.append([(0,)*i + (1,) + (0,)*(n-i-2) + (1,)
                for i in range(1, n-i-1)])
        else:
            rows.append([(0,)+l for l in product([0,1], repeat=n-1)
                if l.count(1) > 1 and l.count(0) < n-2])

    ret = []
    proditer = ([p1,p2,p3] for p1,p2,p3 in product(*rows) if p1!=p2 and p2!=p3 and p1!=p3)
    for M in ProgressBar(maxval=prod(map(len, rows)))(imap(matrix, proditer)):
        if not any(S in ret for S in row_orbit(M)):
            ret.append(M)

    return ret

def checkit(n):
    mats = __get_mats(n-1)
    aF = all_F(mats)
    print sum((i+j+1) * v for (i,j),v in aF.iteritems())
    for k in S.keys():
        print k
        print len(count_structure(n, k))

def Sij(n):
    mats = [(M, M.submatrix(0,0,ncols=n-1)) for M in __get_mats(n)]
    mats_nm1 = __get_mats(n-1)
    ret = []
    for M0,M in mats:
        if any(S in mats_nm1 for S in row_orbit(M)):
                ret.append(M0)
    return ret

