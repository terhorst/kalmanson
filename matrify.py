from functools import wraps
from sage.all import Matrix

def matrify(f):
    @wraps(f)
    def wrapper(*args):
        mat = Matrix(f(*args))
        mat.set_immutable()
        return mat
    return wrapper
