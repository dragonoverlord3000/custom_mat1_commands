import numpy as np
import sympy as sp

def TanKr(V, r, Iu, format="latex", verbose=False):
    """Will calculate the tangential line integral with 
    over the vector field 'V', parameter 'u',
    parameterization 'r(u)' defined in the interval
    u $\in$ [Iu[0], Iu[1]]

    Args:
        V (sp.Matrix): The vector field
        r (sp.Matrix): The parametrization 
        Iu (list or tuple): Interval of parameter
        format (str): The format to return the result in (False -> sympy output)
    """
    x,y,z,u = sp.symbols("x y z u")
    # if rot(V) == 0, then it's a gradient vector field
        
    r_u = sp.diff(r, u)
    
    if len(r) == 2:
        V_r = V.subs({x:r[0], y:r[1]})
    else:
        V_r = V.subs({x:r[0], y:r[1], z:r[2]})
    
    V_r_dot_r_u = (sp.transpose(V_r) * r_u)[0]
    result = sp.integrate(V_r_dot_r_u, (u, Iu[0], Iu[1]))
    
    if verbose:
        pass
    
    if format == "latex":
        return sp.latex(result)
    elif format == "maple":
        return sp.maple_code(result)
    else:
        return result


