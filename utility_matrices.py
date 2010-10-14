import numpy as np
from sage.all import Set, Matrix, vector
from matrify import matrify, symmetrize
   
@symmetrize
@matrify
def V(n, i, j=None):
    i -= 1
    M = np.zeros([n,n], dtype=int)
    if j:
        j -= 1
        M[0:(i+1), (i+1):(j+1)] = 1
        M[(i+1):(j+1), (j+1):(n+1)] = 1
    else:
        M[0:(i+1), (i+1):(n+1)] = 1
    return M

@symmetrize
@matrify
def A(n, i, j=None):
    F = lambda i,j: E(n,i,j)
    if j:
        M = F(i,j+1) + F(i+1,j) - F(i,j) - F(i+1,j+1)
        return M
    else:
        M = F(1,i) + F(i+1,n) - F(1,i+1) - F(i,n)
        return M

    i -= 1
    M = np.zeros([n,n], dtype=int)
    if j:
        j -= 1
        M[0:(i+1), (i+1):(j+1)] = 1
        M[(i+1):(j+1), (j+1):(n+1)] = 1
    else:
        M[0:(i+1), (i+1):(n+1)] = 1
    return M

@symmetrize
@matrify
def E(n, i, j=None):
    M = np.zeros([n,n], dtype=np.int)
    if j:
        M[i-1,j-1] = 1
        return M
    else:
        M[i-1, :] = 1
        M[:, i-1] = 1
        M[i-1,i-1] = 0
        return M

def rays(n):
    s = []
    for i in range(1,n-2):
        for j in range(i+2,n):
            s.append(V(n,i,j))
    for i in range(2,n-1):
        s.append(V(n,i))

    return Set(s)

def lin_space(n):
    return Set([E(n,i+1) for i in range(n)])

def to_Matrix(S):
    n = S[0].nrows()
    inds = np.nonzero(np.triu(np.ones([n,n], dtype=np.int), 1))
    return Matrix(np.vstack([np.array(m)[inds] for m in S]))
