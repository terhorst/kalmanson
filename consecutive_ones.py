from multiprocessing import Pool, cpu_count
from itertools import product, combinations, izip_longest, imap, chain, izip, repeat
from pqtree import PQTree

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def quantify(iterable, pred=bool):
    "Count how many times the predicate is true"
    return sum(imap(pred, iterable))

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def is_circular(tup):
    m,x = tup
    return PQTree(range(m)).ReduceAll(x)

def check_circular_ones(n,m):
    "How many n x m 0/1 matrices possess the circular ones property"
    k = m-1
    rm = range(m)
    p = Pool(cpu_count())
    chunk = 10000
    return sum(p.imap_unordered(is_circular,
            izip(repeat(m),
                combinations((p for p in powerset(range(k))
                    if len(p) > 1 and len(p) < m - 1), n)
                ),
            chunk))
 
def fvector(n):
    for k in range(1,n*(n-3)/2 + 1):
        print check_circular_ones(k,n)
