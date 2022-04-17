import sympy as sp
from tools import cross_product
vector = lambda *args: sp.Matrix([*args])

def vector_potential(W, format="latex", verbose=False):
    """Find a vector-potential of W
    (Dansk: stamvektorfelt) i.e. the vector 
    field V such that rot(W) = V

    note: div(V) = 0 iff rot(W) = V

    Args:
        W (sp.Matrix): Input vector field
    """
    x,y,z,u = sp.symbols("x y z u")
    
    V = -cross_product(vector(x,y,z), sp.integrate(u*W.subs({x:u*x,y:y*u,z:z*u}), (u, 0, 1)))
    if verbose:
        # Is it a vector potential ??? i.e. is div(V) == 0
        pass
    
    if format == "latex":
        return sp.latex(V)
    elif format == "maple":
        return sp.maple_code(V)
    else:
        return V
    


