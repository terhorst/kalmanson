from functools import wraps
from sage.all import Matrix

def matrify(f):
    @wraps(f)
    def wrapper(*args):
        mat = Matrix(f(*args))
        mat.set_immutable()
        return mat
    return wrapper

def sortedmat(f):
    @matrify
    @wraps(f)
    def wrapper(*args):
        return sorted(f(*args))
    return wrapper

def symmetrize(f):
    @wraps(f)
    def wrapper(*args):
        mat = f(*args)
        mat += mat.transpose()
        mat.set_immutable()
        return mat
    return wrapper
