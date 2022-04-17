import sympy as sp

theta = sp.symbols("theta")
Qx = sp.Matrix([[1,0,0],[0, sp.cos(theta), -sp.sin(theta)],[0, sp.sin(theta), sp.cos(theta)]])
Qy = sp.Matrix([[sp.cos(theta), 0, sp.sin(theta)],[0,1,0],[-sp.sin(theta),0,sp.cos(theta)]])
Qz = sp.Matrix([[sp.cos(theta), -sp.sin(theta),0], [sp.sin(theta),sp.cos(theta),0],[0,0,1]])


cross_product = lambda u,v: sp.Matrix([u[1]*v[2] - v[1]*u[2],-u[0]*v[2] + v[0]*u[2],u[0]*v[1] - v[0]*u[1]])

def gradient(F, var1=None, var2=None, var3=None):
    """Returns the gradient of a vector function F
    wrt. x, y and z unless var1 and var2 are specified

    Args:
        v (sp.Matrix): _description_
        var1 (sp.Symbol, optional): _description_. Defaults to None.
        var2 (sp.Symbol, optional): _description_. Defaults to None.
        var3 (sp.Symbol, optional): _description_. Defaults to None.
    """
    x,y,z = sp.symbols("x y z")
    try:    
        if len(F) == 2:
            if var1 and var2:
                return sp.Matrix([sp.diff(F[0], var1), sp.diff(F[1], var2)])
            else:
                return sp.Matrix([sp.diff(F[0], x), sp.diff(F[1], y)])
            
        elif len(F) == 3:
            if var1 and var2 and var3:
                return sp.Matrix([sp.diff(F[0], var1), sp.diff(F[1], var2), sp.diff(F[2], var3)])
            else:
                return sp.Matrix([sp.diff(F[0], x), sp.diff(F[1], y), sp.diff(F[2], z)])
    except Exception as e:
        return sp.Matrix([sp.diff(F,x), sp.diff(F,y), sp.diff(F,z)])
        

def div(V, verbose=False):
    """Returns the divergence of a vector field 'V'
    wrt. 'x', 'y' and 'z'

    Args:
        V (sp.Matrix): _description_
    """
    return sum(gradient(V))
    

def rot(V, verbose=False):
    """Returns the curl of a vector field 'V'
    wrt. 'x', 'y' and 'z'

    Args:
        V (sp.Matrix): _description_
    """
    x,y,z = sp.symbols("x y z")
    return sp.Matrix([
        sp.diff(V[2],y) - sp.diff(V[1],z), 
        -sp.diff(V[2],x) + sp.diff(V[0],z),
        sp.diff(V[1],x) - sp.diff(V[0],y)])

