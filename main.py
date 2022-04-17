# RUN: from main import *
# Then run help(whatever)
# Then run whatever

# Example:
#   >>>from main import *
#   >>>???

# To do:
# 1. Create an exam jupyter notebook
#
# 2.1 Create a numeric plotter for flow curves
# 2.2 Create a numeric plotter for the deformation of lines by a vector field
# 2.3 Create a numeric plotter for the deformation of surfaces by a vector field
# 2.4 Create a numeric plotter for the deformation of solids by a vector field
#
# 3. Create explanations (verbose=True) for all functions
#
# 4. Automatically create parametrization of a surface 
#   4.1 given points in some plane
#   4.2 given description like 'sphere'
#   4.3 given parametrization of solid (then the surface is the rim/edge of the solid)

import sympy as sp
import numpy as np

# Flux = the orthogonal surface integral
from flux_calc import flux
# The curve integral
from Kr_calc import Kr
# The tangential curve integral
from TanKr_calc import TanKr
# The plane integral
from Pr_calc import Pr
# The surface integral
from Fr_calc import Fr
# The volume integral
from Vr_calc import Vr

# Find flow curve from vector field
from flow_curve_calc import flow_curve

# Jacobian calculator
from Jacobian_calc import Jacobian

# Vector field potential calculator
from vector_potential_calc import vector_potential

# Import plotters
from plotters.param_surface_plotter import param_surface_plotter
from plotters.param_curve_plotter import param_curve_plotter

# Show equations nicely using matplotlib
from show_math import show_math

# Random tools
from tools import cross_product, Qx, Qy, Qz, gradient, div, rot

# Sympy helper stuff
Matrix = sp.Matrix
vector = lambda *args: sp.Matrix([*args])
hessian = sp.hessian # e.g. >>>sp.hessian(x**2 + y**2 + 2*z**2 + z, [x,y,z])

# Mathematical functions - sympy
sqrt = sp.sqrt
cos = sp.cos
sin = sp.sin
tan = sp.tan
exp = sp.exp
ln = sp.ln

# Mathematical constants - sympy
pi = sp.pi

# Sympy symbols
u,v,w,x,y,z = sp.symbols("u v w x y z")



















