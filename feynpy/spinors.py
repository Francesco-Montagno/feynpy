import sympy as sp

I = sp.I

def OmegaPlus(En, p):
    return sp.sqrt(En + p)

def OmegaMinus(En, p):
    return sp.sqrt(En - p)

def ChiPlus(theta, phi):
    return sp.Matrix([
        sp.cos(theta/2),
        sp.exp(I*phi)*sp.sin(theta/2)
    ])

def ChiMinus(theta, phi):
    return sp.Matrix([
        -sp.exp(-I*phi)*sp.sin(theta/2),
        sp.cos(theta/2)
    ])

def u_spinor(h, En, p, theta, phi):
    if h == 1/2:
        h = sp.Rational(1,2)
        chi = ChiPlus(theta, phi)
        return sp.Matrix([
            OmegaMinus(En, p)*chi[0],
            OmegaMinus(En, p)*chi[1],
            OmegaPlus(En, p)*chi[0],
            OmegaPlus(En, p)*chi[1],
        ])

    elif h == -1/2:
        h = sp.Rational(-1,2)
        chi = ChiMinus(theta, phi)
        return sp.Matrix([
            OmegaPlus(En, p)*chi[0],
            OmegaPlus(En, p)*chi[1],
            OmegaMinus(En, p)*chi[0],
            OmegaMinus(En, p)*chi[1],
        ])

    else:
        raise ValueError("Helicity must be ±1/2")

def v_spinor(h, En, p, theta, phi):
    if h == 1/2:
        h = sp.Rational(1,2)
        chi = ChiMinus(theta, phi)
        return sp.Matrix([
            OmegaPlus(En, p)*chi[0],
            OmegaPlus(En, p)*chi[1],
            -OmegaMinus(En, p)*chi[0],
            -OmegaMinus(En, p)*chi[1],
        ])

    elif h == -1/2:
        h = sp.Rational(-1,2)
        chi = ChiPlus(theta, phi)
        return sp.Matrix([
            -OmegaMinus(En, p)*chi[0],
            -OmegaMinus(En, p)*chi[1],
            OmegaPlus(En, p)*chi[0],
            OmegaPlus(En, p)*chi[1],
        ])

    else:
        raise ValueError("Helicity must be ±1/2")
