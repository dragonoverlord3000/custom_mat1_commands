# https://matplotlib.org/3.5.1/gallery/mplot3d/lines3d.html

import sympy as sp
import numpy as np

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as Axes3D

def param_curve_plotter(r, Iu, view=None, V=None, show=True):
    """_summary_

    Args:
        r (sp.Matrix): parametric equation in the variable 'u'
        Iu (list): _description_
        view (list, optional): _description_. Defaults to None.
        V (sp.Matrix, optional): _description_. Defaults to None.
        show (bool, optional): _description_. Defaults to True.
    """
    u = sp.symbols("u")
    r_x = sp.lambdify(u, r[0])
    r_y = sp.lambdify(u, r[1])
    r_z = sp.lambdify(u, r[2])

    # Setup data
    u_plot = np.linspace(Iu[0], Iu[1], 400)

    # Setup figure    
    ax = plt.figure().add_subplot(projection='3d')
    
    # Add labels
    ax.set_xlabel("X", fontweight="bold", fontsize=14)
    ax.set_ylabel("Y", fontweight="bold", fontsize=14)
    ax.set_zlabel("Z", fontweight="bold", fontsize=14)

    # Add title
    ax.set_title("Parametric Surface", fontweight="bold", fontsize=16)

    # Adjust figure axes
    if isinstance(view, list):
        ax.set_xlim(view[0][0], view[0][1])
        ax.set_ylim(view[1][0], view[1][1])
        ax.set_zlim(view[2][0], view[2][1])
    
    # Just like maple's 'scaling = "constrained"'
    elif view == "constrained":
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        zlim = ax.get_zlim()
        
        val = max([xlim[1] - xlim[0], ylim[1] - ylim[0], zlim[1] - zlim[0]])
        n_xlim = (val - (xlim[1] - xlim[0]))/2
        ax.set_xlim(xlim[0] - n_xlim, xlim[1] + n_xlim)
        n_ylim = (val - (ylim[1] - ylim[0]))/2
        ax.set_ylim(ylim[0] - n_ylim, ylim[1] + n_ylim)
        n_zlim = (val - (zlim[1] - zlim[0]))/2
        ax.set_zlim(zlim[0] - n_zlim, zlim[1] + n_zlim) 
    
    # Create numerical data
    x = r_x(u_plot)
    y = r_y(u_plot)
    z = r_z(u_plot)

    # Plot data
    ax.plot(x,y,z, label="Parametric Curve")
    ax.legend()
    
    # Show plot
    if show:
        plt.show()


