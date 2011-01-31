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
II = map(vector, [(1,0,0),(0,1,0),(0,0,1),(1,1,1)])

S = {
        'S1': (1, 0, 'X'),
        'S2': (1, 1, 0),
        'S3': (1, 0, 0),
        'S4': (1, 'X', 'X'),
        'S5': (0, 0, 'X'),
        'S6': (1, 1, 1),
        'S7': (0, 0, 0),
        'S8': (0, 'X', 'X'),
        'S9': (1, 1, 'X')
        }

def mat_to_word(ss):
    "".join(map(ALPHA.get, ss.columns()))

def __get_mats(n):
    return [ss_to_matrix(n, ss) for ss in circular_ones(n, 3)]

@memoized
def F(i,j,mats):
    # mats = __get_mats(n)
    return sum( sum(vector(S) in M.columns() for S in I)==i and
                sum(vector(S) in M.columns() for S in II)==j
                for M in mats )

def all_F(mats):
    return dict(((i,j),F(i,j,mats)) for i in range(3) for j in range(4))

def count_structure(n, sty="S1"):
    pat = S[sty]
    rows = []
    for i in pat:
        if i==1:
            rows.append([(0,) + ((1,)*(n-2)) + (0,)])
        elif i==0:
            rows.append([(0,)+l+(1,) for l in permutations((0,)*(n-3) + (1,))])
        else:
            rows.append([(0,)+l for l in product([0,1], repeat=n-1)
                if l.count(1) > 1 and l.count(0) < n-2])

    print rows

    ret = []
    proditer = ([p1,p2,p3] for p1,p2,p3 in product(*rows) if p1!=p2 and p2!=p3 and p1!=p3)
    for M in ProgressBar(maxval=prod(map(len, rows)))(imap(matrix, proditer)):
        if not any(S in ret for S in row_orbit(M)):
            ret.append(M)

    return ret
