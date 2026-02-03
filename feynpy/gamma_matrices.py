import sympy as sp

I = sp.I
Id = sp.eye(4)
_GAMMA_UP = [
    sp.Matrix([
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 0, 0, 0],
        [0, 1, 0, 0]
    ]),
    sp.Matrix([
        [0, 0, 0, 1],
        [0, 0, 1, 0],
        [0, -1, 0, 0],
        [-1, 0, 0, 0]
    ]),
    sp.Matrix([
        [0, 0, 0, -I],
        [0, 0, I, 0],
        [0, I, 0, 0],
        [-I, 0, 0, 0]
    ]),
    sp.Matrix([
        [0, 0, 1, 0],
        [0, 0, 0, -1],
        [-1, 0, 0, 0],
        [0, 1, 0, 0]
    ])
]

def GammaUP(mu: int):
    return _GAMMA_UP[mu]

_GAMMA_DOWN = [
    _GAMMA_UP[0],
    -_GAMMA_UP[1],
    -_GAMMA_UP[2],
    -_GAMMA_UP[3]
]

def GammaDOWN(mu: int):
    return _GAMMA_DOWN[mu]


Gamma5 = sp.Matrix([
    [-1, 0, 0, 0],
    [0, -1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
])

Gamma6 = (Id + Gamma5) / 2
Gamma7 = (Id - Gamma5) / 2


def Slash(p):
    """
    p: 4-vector 
    returns γ^μ p_μ 
    """
    result = sp.zeros(4)
    for mu in range(4):
        result += GammaDOWN(mu) * p[mu]
    return result

def Bar(psi):
    """
    psi: spinor
    returns ψ̄ = ψ† γ^0
    """
    return (psi.H * GammaUP(0))