import sympy as sp

def flow_curve(V, P=[*sp.symbols("x0 y0 z0")], t=sp.Symbol("t")):
    """Calculate flow curve of vector field

    Args:
        V (sp.Matrix): The vector field
        P (list, optional): starting point. Defaults to [*sp.symbols("x0 y0 z0")].
        t (sp.Symbol, optional): parameter. Defaults to sp.Symbol("t")
    """
    x,y,z = sp.symbols("x y z")
    x_f, y_f, z_f = sp.Function("x"), sp.Function("y"), sp.Function("z")
    
    # If in the plane
    if len(V) == 2:
        result = sp.dsolve([sp.diff(x_f(t), t) - V[0].subs({x:x_f(t), y:y_f(t), z:z_f(t)}),
                            sp.diff(y_f(t), t) - V[1].subs({x:x_f(t), y:y_f(t), z:z_f(t)})],
                           ics={x_f(0):P[0], y_f(0):P[1]})
    else:
        result = sp.dsolve([sp.diff(x_f(t), t) - V[0].subs({x:x_f(t), y:y_f(t), z:z_f(t)}),
                            sp.diff(y_f(t), t) - V[1].subs({x:x_f(t), y:y_f(t), z:z_f(t)}),
                            sp.diff(z_f(t), t) - V[2].subs({x:x_f(t), y:y_f(t), z:z_f(t)}),
                            ], ics={x_f(0):P[0], y_f(0):P[1], z_f(0):P[2]})

    return result


