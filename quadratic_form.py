import sympy as sp

form_identifier_table = {}

def quadratic_form(f, format="latex", verbose=False):
    """_summary_

    Args:
        f (sympy expression): expression in x,y and possibly z
        format (str, optional): _description_. Defaults to "latex".
        verbose (bool, optional): _description_. Defaults to False.
    """
    x,y,z = sp.symbols("x y z")
    if z in f.atoms():
        pass
    else:
        H = sp.hessian(f, [x,y,z])
        A = H/2
        pass # ???




