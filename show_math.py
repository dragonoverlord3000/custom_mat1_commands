import sympy as sp
import matplotlib.pyplot as plt

def show_math(formula, fontsize=28, color="black", figsize=(10,2)):
    """Visualize equations with matplotlib

    Args:
        formula (Sympy equation): the formula to show
        fontsize (int): the fontsize of the equation
        color (str): the color of the equation
        figsize (tuple): the dimension of the equation
        
    Example:
        >>>from main import *
        >>>show_math(x*y + x**2 + exp(2*pi)/5 + 8**(2**2))
        
    Limitations:
        Can't show stuff involving matrices
    """
    text_kwargs = dict(ha='center', va='center', fontsize=fontsize, color=color)
    
    plt.figure(figsize=figsize)
    plt.axis("off")
    plt.text(0.5,0.5,"$" + sp.latex(formula) + "$", **text_kwargs)
    plt.show()

