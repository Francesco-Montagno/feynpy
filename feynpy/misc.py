import sympy as sp

g = sp.diag(1, -1, -1, -1)
def metric(mu: int, nu: int):
    return g[mu, nu]

def lower_index(v):
    """
    Lowers an index of a 4x1 vector (column matrix)
    v: UP vector
    returns: DOWN vector
    """
    return g * v  

def raise_index(v):
    """
    Raises an index of a 4x1 vector (column matrix)
    v: DOWN vector
    returns: UP vector
    """
    return g * v


def minkowski_dot(a, b):
    """
    SScalar product a·b with metric (+,-,-,-)
    a, b: 4x1 vectors (column matrices)
    """
    return (a.T * g * b)[0]  