import numpy as np
import sympy as sp

from tools import cross_product

def flux(V, r, Iu, Iv, format="latex", verbose=False):
    """Calculates the flux of the vector field V
    through the surface F_r i.e. the surface parametrized 
    by r(u,v)

    Args:
        V (sp.Matrix): The vector field described using sp.symbols("x y z")
        r (sp.Matrix): A parametrization of the surface
        u (sp.Symbol): One variable
        v (sp.Symbol): The other variable
        Iu (list or tuple): Interval of the u variable
        Iv (list or tuple): Interval of the v variable
    """
    u,v = sp.symbols("u v")
    r_u = sp.diff(r,u)
    r_v = sp.diff(r,v)
    x,y,z = sp.symbols("x y z")
    subs_dict = {x: r[0], y:r[1], z:r[2]}
    crossed = cross_product(r_u, r_v)
    crossed_simplified = sp.simplify(crossed)
    dot_product = (sp.Transpose(V.subs(subs_dict)) * crossed_simplified)[0]
    int_u_eval = sp.integrate(dot_product, (u, Iu[0], Iu[1]))
    result = sp.integrate(int_u_eval, (v, Iv[0], Iv[1]))
    result = sp.simplify(result)
    
    if verbose:
        print(r"r_u^\prime(u,v) = " + sp.latex(r_u))
        print(r"r_v^\prime(u,v) = " + sp.latex(r_v))
        print(r"r_u \times r_v = " + sp.latex(crossed_simplified))
        explanation_str = r"Flux(V, F_r) = \int V \cdot n_{F_r} d\mu = "
        int_v = r"\int_{" + str(Iv[0]) + "}^{" + str(Iv[1]) + "}"
        int_u = r"\int_{" + str(Iu[0]) + "}^{" + str(Iu[1]) + "}"
        subs_dict = {x: r[0], y:r[1], z:r[2]}
        explanation_str += int_v + int_u + sp.latex(V.subs(subs_dict)) + r"\cdot" +r"\left(" + sp.latex(r_u) + r"\times" + sp.latex(r_v) + r"\right)" + "dudv = "
        explanation_str += int_v + int_u + sp.latex(V.subs(subs_dict)) + r"\cdot" + sp.latex(crossed_simplified) + "dudv = "
        explanation_str += int_v + int_u + sp.latex(dot_product) + "dudv = "
        explanation_str += int_v + sp.latex(int_u_eval) + "dv = "
        explanation_str += sp.latex(result)
        print("\n")
        print(explanation_str)
        print("\n")
    
    if format == "latex":
        return sp.latex(result)
    elif format == "maple":
        return sp.maple_code(result)
    else:
        return result
    
    








