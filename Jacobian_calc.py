import sympy as sp
from tools import cross_product


def Jacobian(r, format="latex", verbose=False):
    """Find Jacobian given parameterization

    Args:
        r (sp.Matrix): parametrization of 'u' or 'u' and 'v' or 'u', 'v' and 'w'
        format (str, optional): _description_. Defaults to "latex".
        verbose (bool, optional): _description_. Defaults to False.
    """
    u,v,w = sp.symbols("u v w")
    r_atoms = r.atoms()
    if u in r_atoms:
        if v in r_atoms:
            if w in r_atoms:
                Jacobian = JacobianSolid(r, format=format, verbose=verbose)
            else:
                if len(r) == 2:
                    Jacobian = JacobianPlanar(r, format=format, verbose=verbose)
                else:
                    Jacobian = JacobianSurface(r, format=format, verbose=verbose)
        else:
            Jacobian = JacobianCurve(r, format=format, verbose=verbose)

    Jacobian = sp.simplify(Jacobian)
    if format == "latex":
        return sp.latex(Jacobian)
    elif format == "maple":
        return sp.maple_code(Jacobian)
    else:
        return Jacobian


# Jacobian helper functions
def JacobianSolid(r, format="latex", verbose=False):
    """Jacobian for solid figure

    Args:
        r (_type_): _description_
        format (str, optional): _description_. Defaults to "latex".
        verbose (bool, optional): _description_. Defaults to False.
    """
    u,v,w = sp.symbols("u v w")
    r_u = sp.diff(r, u)
    r_v = sp.diff(r, v)
    r_w = sp.diff(r, w)
    result = sp.Abs(sp.det(sp.Matrix([[r_u[0],r_u[1],r_u[2]], [r_v[0], r_v[1], r_v[2]], [r_w[0], r_w[1], r_w[2]]])))
    
    if verbose:
        print(f"""
              r_u^\\prime(u,v) &= {sp.latex(r_u)} \\\\
              r_v^\\prime(u,v) &= {sp.latex(r_v)} \\\\
              r_w^\\prime(u,v) &= {sp.latex(r_v)} \\\\
              Jakobi_r = \\left| \det\\begin{bmatrix} {r_u[0]} & {r_v[0]} & {r_w[0]} \\\\ {r_u[1]} & {r_v[1]} & {r_w[1]} \\\\ {r_u[2]} & {r_v[2]} & {r_w[2]} \\end{bmatrix} \\right| = {result}
              """)    

    return result

def JacobianSurface(r, format="latex", verbose=False):
    """Jacobian for surface

    Args:
        r (_type_): _description_
        format (str, optional): _description_. Defaults to "latex".
        verbose (bool, optional): _description_. Defaults to False.
    """
    u,v = sp.symbols("u v")
    r_u = sp.diff(r, u)
    r_v = sp.diff(r, v)
    result = sp.Abs(cross_product(r_u, r_v).norm())

    if verbose:
        print(f"""
              r_u^\\prime(u,v) &= {sp.latex(r_u)} \\\\
              r_v^\\prime(u,v) &= {sp.latex(r_v)} \\\\
              r_w^\\prime(u,v) &= {sp.latex(r_v)} \\\\
              Jakobi_r = \\left| \det\\begin{"\{bmatrix\}"} i & j & k \\\\ {r_u[0]} & {r_u[1]} & {r_u[2]} \\\\ {r_v[0]} & {r_v[1]} & {r_v[2]} \\end{"\{bmatrix\}"} \\right| = {sp.latex(result)}
              """)    
    
    return result

def JacobianPlanar(r, verbose=False):
    """Jacobian for planar surface

    Args:
        r (_type_): _description_
        format (str, optional): _description_. Defaults to "latex".
        verbose (bool, optional): _description_. Defaults to False.
    """
    u,v = sp.symbols("u v")
    r_u = sp.diff(r, u)
    r_v = sp.diff(r, v)
    result = sp.Abs(sp.det(sp.Matrix([[r_u[0], r_u[1]], [r_v[0], r_v[1]]])))

    if verbose:
        print(f"""
              r_u^\\prime(u,v) &= {sp.latex(r_u)} \\\\
              r_v^\\prime(u,v) &= {sp.latex(r_v)} \\\\
              Jakobi_r = \\left| \det\\begin{"\{bmatrix\}"} {r_u[0]} & {r_v[0]} \\\\ {r_u[1]} & {r_v[1]} \\end{"\{bmatrix\}"} \\right| = {sp.latex(result)}
              """)
    
    return result


def JacobianCurve(r, format="latex", verbose=False):
    """Jacobian for parametric curve

    Args:
        r (_type_): _description_
        format (str, optional): _description_. Defaults to "latex".
        verbose (bool, optional): _description_. Defaults to False.
    """
    u = sp.Symbol("u")
    diff_u = sp.diff(r,u)
    result = diff_u.norm()
    
    return result



