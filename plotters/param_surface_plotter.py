import sympy as sp
import numpy as np

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as Axes3D

from tools import cross_product


# Inspiration from: https://www.youtube.com/watch?v=l9kAJd4TkaA
def param_surface_plotter(r, Iu, Iv, view=None, V=None, n_vec=None, show=True):
    """_summary_

    Args:
        r (sp.Matrix): parametric equation in 'u' and 'v'
        Iu (list or tuple): _description_
        Iv (list or tuple): _description_
        view (list, optional): set x,y- and z-limits e.g. view=[[0,1],[0,6.28],[2.09, 10]] or view="constrained". Defaults to None.
        V (sp.Matrix, optional): optional vector field to plot. Defaults to None.
        n_vec (list, optional): optional unit normal vector to plot e.g. [u_0,v_0]. Defaults to None.
        show (bool, optional): _description_. Defaults to True.        
    """
    u,v = sp.symbols("u v")
    r_x = sp.lambdify([u,v], r[0])
    r_y = sp.lambdify([u,v], r[1])
    r_z = sp.lambdify([u,v], r[2])

    u_plot, v_plot = np.linspace(Iu[0], Iu[1], 400), np.linspace(Iv[0], Iv[1], 400)
    u_plot, v_plot = np.meshgrid(u_plot, v_plot)

    # Create numerical data
    x = r_x(u_plot, v_plot)
    y = r_y(u_plot, v_plot)
    z = r_z(u_plot, v_plot)

    # Create figure
    fig = plt.figure("Parametric Surfaces")
    ax = fig.add_subplot(111, projection="3d")

    # Plot the data
    h = ax.plot_surface(x, y, z, cmap="jet", edgecolor="k") # black edgecolor

    # Colorbar
    fig.colorbar(h)
    
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
           
    # If vector field: https://matplotlib.org/3.1.0/gallery/mplot3d/quiver3d.html ???
    if n_vec:
        x_0,y_0,z_0 = tuple(r.subs({u:n_vec[0], v:n_vec[1]}))
        r_cross = cross_product(sp.diff(r,u), sp.diff(r,v))
        x_temp, y_temp, z_temp = tuple(r_cross.subs({u:n_vec[0], v:n_vec[1]}))
        length = (x_temp**2 + y_temp**2 + z_temp**2)**0.5
        x_temp, y_temp, z_temp = x_temp/length, y_temp/length, z_temp/length
        x_1, y_1, z_1 = x_0 + x_temp, y_0 + y_temp, z_0 + z_temp

        ax.quiver(x_0,y_0,z_0, x_1, y_1, z_1, color="red", linewidths=2, normalize=True)
    
    # Show the figure
    if show:
        plt.show()




