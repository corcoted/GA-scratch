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
A, B, C, D = symbols("A B C D", real=True)

# %%
M = pga2.lt([[A, C, 0],[B,D,0],[0,0,1]]) # generic ray transfer matrix (note the list is the transpose of the matrix)


# %%
# helper functions
def point(x, y):
    return x*e20 + y*e01 + e12

def line(a, b, c):
    return a*e1+b*e2+c*e0


# %%
# define basic centered rtms

# %%
f, d, R = symbols("f d R", real=True) # focal length, distance, radius
n1 ,n2 = symbols("n1 n2", positive=True) # indices of refraction

# %%
a, b, c = symbols("a b c", real=True)

# %%
M_lens = lambda f: pga2.lt([[1, -1/f, 0],[0,1,0],[0,0,1]])
M_flatref = lambda n1, n2: pga2.lt([[1, 0, 0],[0,n1/n2,0],[0,0,1]])
M_flatmirror = pga2.lt([[1, 0, 0],[0,-1,0],[0,0,1]])
M_sphereref = lambda n1, n2, R: pga2.lt([[1, (n1-n2)/(R*n2), 0],[0,n1/n2,0],[0,0,1]])
M_spheremirror = lambda R: pga2.lt([[1,- 2/R, 0],[0,1,0],[0,0,1]])
M_distance = lambda d: pga2.lt([[1,0,0],[d,1,0],[0,0,1]])


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
# try composing some rtms together
(M_distance(2*f)*M_lens(f)*M_distance(2*f))(point(-x,y))

# %%
# define translation and rotation linear operators
M_trans = lambda x,y: pga2.lt(1+y/2*e20-x/2*e01)
M_rot = lambda a, p=e12: pga2.lt(cos(a/2)-sin(a/2)*p) # rotation by angle a about point p

# %%
M_rot(a)(point(x,y)).trigsimp()

# %%
M_rot(pi,point(0,-c/b))(line(a,b,c))

# %%
simplify(trigsimp(_95.obj))

# %% [markdown]
# One gap to close is to match our work with the 3x3 matrix (Siegman or Tovar) and the 4x4 matrix of Shaomin.

# %%
M_rot(-f)(M(M_rot(f)(line(a,-1,b))))

# %%
series(_56.get_coefs(1)[0],f,n=2)

# %%
series(_56.get_coefs(1)[1],f,n=2)

# %%
series(_56.get_coefs(1)[2],f,n=2)

# %%
series(-_57/_59,f,n=2).expand().collect([a,b])

# %%
series(-_58/_59,f,n=2).expand().collect([a,b])

# %%
series(-_59/_59,f,n=2).expand().collect([a,b])

# %%
# the above confirms that this agrees with Shaomin up to normalization and linear dependence on the ray parameters.

# %%
M_trans(0,-d)(M(M_trans(0,d)(line(a,-1,b))))

# %%
# This agrees with Shaomin eq. 6 for shifted ABCD matrix!

# %% [markdown]
# Tilt-shift example
#
# Show that the object plane, image plane, and focal plane of a thin lens meet at at point

# %%
# start with vertical object plane
obj=line(1,0,-d)

# %%
obj

# %%
img=M_lens(f)(obj)

# %%
obj^img # their intersection point is infinitely far away in the y direction.


# %%
def J(x):
    # J map for pga2 multivectors
    coef_list = x.blade_coefs()
    # flip sign of e02 and e1
    coef_list[2] = -coef_list[2]
    coef_list[5] = -coef_list[5]
    return pga2.mv(sum(coef_list[7-mm]*flatten(pga2.blades)[mm] for mm in range(8)))


# %%
J(J(point(x,y)))

# %%
J(J(line(a,b,c)))


# %%
def vee(A, B):
    # regressive product
    return J(J(A)^J(B))


# %%
#now tilt the object plane about where it crosses the axis
obj2=vee(point(-d,0),(sin(a)*e20+cos(a)*e01))

# %%
img2 = M_lens(f)(obj2)

# %%
img2

# %%
obj2^img2

# %%
lensplane = line(1,0,0) # y axis

# %%
lensplane

# %%
obj2^img2^lensplane == 0 # do these all intersect at the same point?

# %%
Euc_normalize(obj2^img2)

# %%
# This is the point of intersection, x=0, y=the first coefficient

# %%
img2^e2

# %%
Euc_normalize(img2^e2)

# %%
simplify((1/f-1/d)**-1)
