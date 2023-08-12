# ---
# title: 2D PGA in `galgebra`
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-autoscroll,-collapsed,-scrolled,-trusted,-ExecuteTime
#     notebook_metadata_filter: kernelspec,jupytext,-jupytext.text_representation.jupytext_version
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python [conda env:galgebra]
#     language: python
#     name: conda-env-galgebra-py
# ---

# %% [markdown]
# Here I'm trying to do some symbolic calculations in 2D PGA using the `galgebra` package.

# %%
from sympy import *
from galgebra.ga import Ga
from galgebra.printer import latex
from IPython.display import Math

# tell sympy to use our printing by default
init_printing()

# %% [markdown]
# Attempt to generate the algebra for PGA directly:

# %%
wxy = (w, x, y) = symbols('w x y', real=True)
pga2 = Ga('e_0 e_1 e_2', g=[0,1, 1], coords=wxy)
grad = pga2.grad

# %%
# not sure of the best choice here, see docs: 
# https://galgebra.readthedocs.io/en/latest/generated/galgebra.ga.html?galgebra.ga.Ga.dual_mode
pga2.dual_mode(mode='+I')

# %%
e0, e1, e2 = pga2.mv()

# %%
# define some aliases to simplify writing
e01 = e0*e1; e20=e2*e0; e12=e1*e2; e012=e0*e1*e2


# %%
# helper functions
def point(x, y):
    return x*e20 + y*e01 + e12

def line(a, b, c):
    return a*e1+b*e2+c*e0


# %%
# Euclidean norm
def Euc_norm(p):
    if p.pure_grade() < 0: #not single grade
        return nan
    elif p.pure_grade() == 0:
        return p.get_coefs(0)
    elif p.pure_grade() == 3:
        return 0
    elif p.pure_grade() == 1: # line
        return sqrt(p.get_coefs(1)[1]**2 + p.get_coefs(1)[2]**2)
    elif p.pure_grade() == 2: # point
        return p.get_coefs(2)[2]


# %%
# Ideal norm
def Ideal_norm(p):
    if p.pure_grade() < 0: #not single grade
        return nan
    elif p.pure_grade() == 0:
        return 0
    elif p.pure_grade() == 3:
        return p.get_coefs(3)
    elif p.pure_grade() == 1: # line
        return p.get_coefs(1)[0]
    elif p.pure_grade() == 2: # point
        return sqrt(p.get_coefs(2)[0]**2 + p.get_coefs(2)[1]**2)


# %%
# Euc_normalize
def Euc_normalize(p):
    n = Euc_norm(p)
    if n==nan:
        return nan
    return p/n


# %%
# define translation and rotation linear operators
M_trans = lambda x,y: pga2.lt(1+y/2*e20-x/2*e01)
M_rot = lambda a, p=e12: pga2.lt(cos(a/2)-sin(a/2)*p) # rotation by angle a about point p

# %%
L, px, py = symbols("L, px, py", real=True)

# %%
s = e0 + x*e1+y*e2+L*e12 + px*e20 + py*e01

# %%
s

# %%
dx, dy, ω = symbols("dx, dy, ω", real=True)

# %%
(1+dx/2*e20 + dy/2*e01)*s*(1-dx/2*e20 -dy/2*e01)

# %%
sd = L*e0 + px*e1+py*e2+e12 + x*e20 + y*e01

# %%
(1+dy/2*e20 - dx/2*e01)*sd*(1-dy/2*e20 +dx/2*e01)

# %%
from galgebra.mv import Mv

# %%
out = Mv(ga=pga2)
print([i for i in s.blade_coefs()])

# %%
(~s).blade_coefs()[::-1]

# %%
pga2.indexes_lst
