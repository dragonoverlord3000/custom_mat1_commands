import numpy as np
import sympy as sp

def Vr(r, Iu, Iv, Iw, f = lambda x,y,z: 1, format="latex", verbose=False):
    """_summary_

    Args:
        r (_type_): _description_
        Iu (_type_): _description_
        Iv (_type_): _description_
        Iw (_type_): _description_
        f (function, optional): density function. Defaults to lambdax.
        format (str, optional): _description_. Defaults to "latex".
        verbose (bool, optional): _description_. Defaults to False.
    """
    x,y,z,u,v,w = sp.symbols("x y z u v w")
    r_u = sp.diff(r,u)
    r_v = sp.diff(r,v)
    r_w = sp.diff(r,w)
    Jacobian = sp.det(sp.Matrix([[r_u[0], r_u[1], r_u[2]], [r_v[0], r_v[1], r_v[2]], [r_w[0], r_w[1], r_w[2]]]))
    int_u = sp.integrate(f(r[0],r[1],r[2]) * Jacobian, (u,Iu[0],Iu[1]))
    int_v = sp.integrate(int_u, (v,Iv[0],Iv[1]))
    result = sp.integrate(int_v, (w, Iw[0], Iw[1]))

    if verbose:
        pass
    
    return result



