import sympy as sp

En, p, θ, φ = sp.symbols('En p theta phi', real=True)

def MomUP(En, p, θ, φ):
    return sp.Matrix([
        En,
        p*sp.sin(θ)*sp.cos(φ),
        p*sp.sin(θ)*sp.sin(φ),
        p*sp.cos(θ)
    ])

def MomDOWN(En, p, θ, φ):
    return sp.Matrix([
        En,
        -p*sp.sin(θ)*sp.cos(φ),
        -p*sp.sin(θ)*sp.sin(φ),
        -p*sp.cos(θ)
    ])
