import sympy as sp
from tools import rot
from TanKr_calc import TanKr

def Cirk(V, r, Iu, format="latex", verbose=False):
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
    # if rot(V) == 0, then it's a gradient vector field
    if not any(rot(V, verbose=verbose)):
        if verbose:
            pass
        result = 0
    # elif stokes is usefull
    else:
        result = TanKr(V, r, Iu, format=format, verbose=verbose)
    
    if format == "latex":
        return sp.latex(result)
    elif format == "maple":
        return sp.maple_code(result)
    else:
        return result


