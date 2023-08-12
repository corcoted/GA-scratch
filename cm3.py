# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-autoscroll,-collapsed,-scrolled,-trusted,-ExecuteTime
#     notebook_metadata_filter: kernelspec,jupytext,-jupytext.text_representation.jupytext_version
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python (galgebra)
#     language: python
#     name: galgebra
# ---

# %% [markdown]
# <h2>*Introduction*</h2>
#
# This file enables a user to construct and manipulate geometric objects in $\mathbb{R}^3$. The constructions and manipulations are performed using a conformal model of $\mathbb{R}^3$. A user need not know much about the conformal model, as all constructions and manipulations are via the functions provided here.   
#
# My intent is that the functions are self documenting through their code and comments. 
#
# Vectors passed to the functions must be in conformal representation. 
# Exceptions: pt, which converts a 3D point to a conformal point, a 3D normal vector $\mathbf{n}$, and a 3D parallel bivector $\mathbf{B}$.
#
# Objects returned by the functions are in conformal representation. Exception: tp, which returns a 3D point.
#
# To start, pull down the "Cell" menu item and choose "Run All".
#
# Comments and proposed changes or additions are welcome.

# %%
# Conformal Model, Amsterdam convention.  Dorst et al. p. 361
from sympy import *
from galgebra.ga import Ga
from galgebra.mv import *
# from lt import *
# from sympy import *

cm3coords = (o,x,y,z,infty) = symbols('o 1 2 3 infty', real=True)
cm3g = '0 0 0 0 -1, 0 1 0 0 0, 0 0 1 0 0, 0 0 0 1 0, -1 0 0 0 0'
cm3 = Ga('o e_1 e_2 e_3 oo', g = cm3g, coords = cm3coords)
(eo, e1, e2, e3, eoo) = cm3.mv()
ep = eo - eoo/2  # ep^2 = +1  GACS 408
em = eo + eoo/2  # em^2 = -1
E = eo^eoo
Ga.dual_mode('Iinv+')
#cm3coords = (o,x,y,z,infty) = symbols('o x y z \infty', real=True)
#cm3 = Ga('o e_x e_y e_z \infty', g = cf3g, coords = cf3coords)
from IPython.display import display


# %%
def pt(arg): # R^3 vector --> conformal point. 
    if isinstance(arg,str):           # Return general 3D point
        v = cm3.mv(arg, 'vector')     # General conformal vector 
        v = v + (v < eoo)*eo + (v < eo)*eoo  # 3D part 
        v = eo + v + (v<v)*eoo/2
    elif arg == 0:
        v = eo
    elif (arg < eoo) == 0:    # Return point for 3D vector in arg
        v = eo + arg + (arg<arg)*eoo/2
    else: v = arg     # arg already in conformal representation   
    return(v)


# %%
def tp(arg): # conformal point --> R^3 vector
    if isinstance(arg,str):   # Return general 3D vector
        v = cm3.mv(arg, 'vector')
    else:                     # Return 3D vector part of arg
        v = arg
    v = v + (v < eoo)*eo + (v < eo)*eoo
    return(v)       


# %%
def normalize(v): 
    if (v < eoo) == 0: # Normalize 3D vector
        return(v/sqrt((v<v).scalar()))
    else:  # Normalize conformal vector: set eo coeff to 1.
        return(-v/(v<eoo))


# %%
def scalar(arg):
    return(cm3.mv(arg, 'scalar')) # Save user from typing all this


# %% [markdown]
# <h4>* Create direct representations of geometric objects *</h4>

# %%
def round(*args):  # args are conformal points
    ans = args[0]
    for i in range(1,len(args)):
        ans = ans ^ args[i]
    return(ans)


# %%
def flat(*args):   # args are conformal points
    return(round(*args) ^ eoo)


# %%
def line(p,q): # If q is 3D, line thru p parallel to q returned
    return(flat(p,q))


# %%
def plane(p,q,r):
    return(flat(p,q,r))


# %%
def circle(p,q,r):
    return(round(p,q,r))


# %%
def sphere(p,q,r,s):
    return(round(p,q,r,s))


# %% [markdown]
# <h4>* Create dual representations of geometric objects *</h4>

# %%
def dualLine(p, B): # thru point p, orthogonal to 3D bivector B
    return(p < (B*eoo))  # A vector


# %%
def dualPlane(p,n):    # n: GA^3 normal vector    
    m = normalize(n)
    if isinstance(p,(int, long, float)):
        p = scalar(p)         # Python scalar -> GAlgebra scalar
    if (p!=0) and ((p<p)==0): # p: point on plane. 
        return(p < (m^eoo))   # a vector
    else:                     # p: distance to origin.
        return(m + (p*eoo))   # a vector


# %%
def dualSphere(c,rho):  # c:center. 
    if isinstance(rho,(int, long, float)):
        rho = scalar(rho)   # Python scalar -> GAlgebra scalar
    if (rho!=0) and ((rho<rho)==0):  # rho: point on sphere 
        return(rho < (c ^ eoo))  
    else:                            # rho: radius. 
        return(c - (rho*rho*eoo)/2)  # A vector     


# %%
def dualCircle(c,rho,n): # c:center. rho:radius. n:normal vector
    ds = dualSphere(c,rho)
    dp = dualPlane(c,n)    
    return(ds^dp)          # A BIvector 


# %% [markdown]
# <h4>*  Geometric operations *</h4>

# %%
def translate(object,a3): # a3: 3D vector
    return(1 - a3*eoo/2)*object*(1 + a3*eoo/2)


# %%
def rotate(object,itheta):
    return(exp(-itheta/2)*object*exp(itheta/2))


# %%
def invert(p, norm=False):   # GACS 513
    ans = -(eo - eoo/2)*p*(eo - eoo/2) 
    if norm:
        ans = normalize(ans)
    return(ans)


# %%
# Reflect point p in hyperplane with normal 3D vector n.
def reflect(p,n):
    return(-n*p*(n/norm2(n)))  


# %%
# Can be considerably simplified: A Covariant Approach ..., 16 
def dilate(p, alpha, norm = False):  # Dilate by alpha (> 0)
    ans = exp(E*ln(alpha)/2)*p*exp(-E*ln(alpha)/2)
    if norm:
        ans = normalize(ans)
    return(ans)


# %% [markdown]
# <h4>* Play *</h4>

# %%
pt(0*e1)

# %%
translate(_,1*e1)

# %%
E

# %%
translate(1*E,1*e1)

# %%
translate(eo*e1*e2*e3*eoo,1*e1)

# %%
f = symbols("f", real=True)
