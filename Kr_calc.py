import numpy as np
import sympy as sp

def Kr(r,Iu, f= lambda x,y,z: 1, format="latex", verbose=False):
    """Will calculate the line integral with 
    density function 'f(x,y,z)', parameter 'u',
    parameterization 'r(u)' defined in the interval
    u $\in$ [Iu[0], Iu[1]]

    Args:
        r (sp.Matrix): The parametrization
        Iu (list or tuple): Interval of parameter
        f (function): The density function dependent on 3 variables
        format (str): The format to return the result in (False -> sympy output)
    """
    u = sp.Symbol("u")
    r_u = sp.diff(r, u)
    Jacobian = r_u.norm()
    if len(r) == 2:
        result = sp.integrate(f(r[0],r[1]) * Jacobian,(u, Iu[0], Iu[1]))
    else:
        result = sp.integrate(f(r[0],r[1],r[2]) * Jacobian,(u, Iu[0], Iu[1]))
        
    
    if verbose:
        print(r"r_u^\prime = " + sp.latex(r_u))
        print(r"\textrm{Jacobian}_r = " + sp.latex(Jacobian))
        print(r"\int_{" + str(Iu[0]) + r"}^{" + str(Iu[1]) + "}" + sp.latex(f(r[0],r[1],r[2]) * Jacobian) + "du = " + sp.latex(result))
    
    if format == "latex":
        return sp.latex(result)
    elif format == "maple":
        return sp.maple_code(result)
    else:
        return result


