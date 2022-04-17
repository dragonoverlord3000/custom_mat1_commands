import numpy as np
import sympy as sp

def Pr(r, Iu, Iv, f = lambda x,y: 1, format="latex", verbose=False):
    """Will calculate the planer integral (Dansk: plan integrale)
    of 'r' wrt. 'u' and 'v' over the intervals 'Iu' and 'Iv' with
    density function 'f' using a format specified by 'format'
    
    - note that this is just a special case of the surface integral

    Args:
        r (sp.Matrix): the parameterization of 'r'
        Iu (list or tuple): the u-interval
        Iv (list or tuple): the v-interval
        f (function, optional): density function. Defaults to lambdax.
        format (str, optional): the output format. Defaults to "latex".
        verbose (bool, optional): how much to print. Defaults to False.
    """
    u,v,x,y = sp.symbols("u v x y")
    r_u = sp.diff(r, u)
    r_v = sp.diff(r, v)
    Jacobian = sp.det(sp.Matrix([[r_u[0], r_u[1]], [r_v[0], r_v[1]]]))
    int_u = sp.integrate(f(r[0], r[1]) * Jacobian, (u, Iu[0], Iu[1]))
    result = sp.integrate(int_u, (v, Iv[0], Iv[1]))
    
    if verbose:
        pass
    
    if format == "latex":
        return sp.latex(result)
    elif format == "maple":
        return sp.maple_code(result)
    else:
        return result





