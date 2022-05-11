vector = lambda *args: sp.Matrix([*args])
import sympy as sp


form_identifier_table = {}

def quadratic_form(f, format="latex", verbose=True):
    """_summary_

    Args:
        f (sympy expression): expression in x,y and possibly z with less than or equal to second degree terms
            e.g. x**2 + x*y + 4*y**2
        format (str, optional): _description_. Defaults to "latex".
        verbose (bool, optional): _description_. Defaults to True.
            The explanation is pretty usefull to copy-paste
            
    Returns:
        matrices: 'Q','A','Lambda' such that Q is orthogonal and Q.T * A * Q = Lambda
        equation: 'new_eq' - the equation in the new variables induced by 'Q'
    """
    x,y,z = sp.symbols("x y z")
    if z in f.atoms():
        og_var_vec = vector(x,y,z)
        
        # Setup matrix and diagonalize
        H = sp.hessian(f, [x,y,z])
        A = H/2
        Q, D = A.diagonalize()
        
        # Normalize
        Q1,Q2,Q3 = sp.GramSchmidt([Q.col(0), Q.col(1), Q.col(2)], True)
        Q[:,0] = sp.simplify(Q1)
        Q[:,1] = sp.simplify(Q2)
        Q[:,2] = sp.simplify(Q3)
        Lambda = Q.T * A * Q
        
        # Get linear terms
        linear_x = sp.diff(f,x).subs({x:0,y:0})
        linear_y = sp.diff(f,y).subs({x:0,y:0})
        linear_z = sp.diff(f,z).subs({x:0,y:0})
        linear_vec = vector(linear_x, linear_y, linear_z)
        
        const_term = f.subs({x:0,y:0,z:0})
        
        # Equation in new variables
        x1,y1,z1 = sp.symbols("x1 y1 z1")
        var_vec = vector(x1,y1,z1)
        new_eq = sp.simplify((var_vec.T * Lambda * var_vec)[0] + (var_vec.T * Q * linear_vec)[0])
        
    else:
        og_var_vec = vector(x,y)
        
        # Setup matrix and diagonalize
        H = sp.hessian(f, [x,y])
        A = H/2
        Q, D = A.diagonalize()
        
        # Normalize
        Q1,Q2 = sp.GramSchmidt([Q.col(0), Q.col(1)], True)
        Q[:,0] = Q1
        Q[:,1] = Q2
        Lambda = Q.T * A * Q
        
        # Get linear terms
        linear_x = sp.diff(f,x).subs({x:0,y:0})
        linear_y = sp.diff(f,y).subs({x:0,y:0})
        linear_vec = vector(linear_x, linear_y)

        const_term = f.subs({x:0,y:0})
        
        # Equation in new variables
        x1,y1 = sp.symbols("x1 y1")
        var_vec = vector(x1,y1)
        new_eq = sp.simplify((var_vec.T * Lambda * var_vec)[0] + (var_vec.T * Q * linear_vec)[0])
        
    if verbose:
        print("\n\n--------------- Derivation ---------------")
        print(f"""
              f(x,y{"" if len(og_var_vec) <= 2 else ",z"}) &= {sp.latex(f)} = {sp.latex(og_var_vec.T)}{sp.latex(A)}{sp.latex(og_var_vec)} + {sp.latex(og_var_vec.T)}{sp.latex(linear_vec)} + {const_term} \\\\ 
              &= {sp.latex(var_vec.T)}{sp.latex(Lambda)}{sp.latex(var_vec)} + {sp.latex(var_vec.T)}{sp.latex(Q * linear_vec)} + {const_term} \\\\
              &= {sp.latex(new_eq)}
              """)
        print("---------------------------------------------")
        
    return (Q, A, Lambda, new_eq)   



