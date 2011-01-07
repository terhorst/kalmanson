from multiprocessing import Pool, cpu_count
import operator
from itertools import *
from pqtree import PQTree
from math import factorial
from memoized import memoized

def binomial(n,k):
    return factorial(n) / factorial(k) / factorial(n-k)

def mypowerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(2, len(s)))

def reduce_all(tup):
    m,x = tup
    x = tuple(x)
    return (x, PQTree(range(m)).ReduceAll(x))

def check_consecutive_ones(n,m):
    rm = range(m)
    p = Pool(cpu_count())
    chunk = 10000
    return sum(p.imap_unordered(reduce_all,
            izip(repeat(m), combinations((p for p in powerset(range(m))), n)),
            chunk))

@memoized
def circular_ones(n,m,non=False):
    "Which nXm 0/1 matrices possess the circular ones property"
    p = Pool(cpu_count())
    chunk = 10000

    if False and m>1 and (num_circular_ones(n,m-1)*num_circular_ones(n,1) <
            binomial(num_circular_ones(n,1), m)):
        poss = tuple(set(ifilter(lambda s: len(s) == m,
                starmap(operator.add, product(circular_ones(n,m-1), circular_ones(n,1))))))
    else:
        poss = combinations(mypowerset(range(n-1)), m)

    return tuple(set([ss for ss,circ in p.imap(reduce_all, izip(repeat(n), poss))
            if (not circ if non else circ)]))

def set_union(a,b):
    return a.union(b)

@memoized
def num_circular_ones(m,n):
    return len(circular_ones(m,n))

def fvector(n):
    for k in range(1, n*(n-3)/2 + 1):
        print num_circular_ones(n,k)
