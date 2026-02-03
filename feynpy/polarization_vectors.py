import sympy as sp

I = sp.I

def EpsilonUP(lam, En, p, theta, phi):
    if lam == 1:
        return (sp.exp(I*phi)/sp.sqrt(2)) * sp.Matrix([
            0,
            -sp.cos(theta)*sp.cos(phi) + I*sp.sin(phi),
            -sp.cos(theta)*sp.sin(phi) - I*sp.cos(phi),
            sp.sin(theta)
        ])

    elif lam == -1:
        return -(sp.exp(-I*phi)/sp.sqrt(2)) * sp.Matrix([
            0,
            -sp.cos(theta)*sp.cos(phi) - I*sp.sin(phi),
            -sp.cos(theta)*sp.sin(phi) + I*sp.cos(phi),
            sp.sin(theta)
        ])

    elif lam == 0:
        return (1/sp.sqrt(En**2 - p**2)) * sp.Matrix([
            p,
            En*sp.sin(theta)*sp.cos(phi),
            En*sp.sin(theta)*sp.sin(phi),
            En*sp.cos(theta)
        ])

    else:
        raise ValueError("Helicity must be -1, 0, or +1")
