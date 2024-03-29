# ---
# title: Brief introduction to linear algebra using PGA
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
# Ted Corcovilos, 2021-03-18

# %%
# import libraries
from sympy import *
from galgebra.ga import Ga
# setup notebook pretty printing
init_printing()

# %%
# set up the algebra
pga3coords = (w,x,y,z) = symbols('w x y z', real=True)
pga3 = Ga('e_0 e_1 e_2 e_3',g=[0,1,1,1], coords=pga3coords)

e0, e1, e2, e3 = pga3.mv()


# %%
# define some useful functions for PGA
def J3(x):
    # J map for pga3 multivectors to calculate the orthogonal complement
    coef_list = x.blade_coefs()
    # flip sign of e02 and e1
    signs = [1,
            1,-1,1,-1,
            1,-1,1,1,-1,1,
            1,-1,1,-1,
            1]
    size = len(signs)
    return pga3.mv(sum(signs[mm]*coef_list[mm]*flatten(pga3.blades)[size-1-mm] for mm in range(size)))

def vee3(a,b):
    # regressive product for pga3
    return J3(J3(a)^J3(b))

def norm(x):
    # norm.  Note that this does not keep the sign.
    return sqrt(((~x)*x).scalar())

def norm_ideal(x):
    # ideal norm.  Note that this does not keep the sign.
    return sqrt((~J3(x)*J3(x)).scalar())

def normalize(x):
    # normalize a MV
    return x/norm(x)


# %%
t= symbols('t', real=True) # real parameter to use later

# %% [markdown]
# ## Systems of linear equations
# Each linear equation may be represented by a vector in PGA.
#
# Let's solve the system
# $$
# \begin{aligned}
# a:&& x + \phantom{2}y + z -1 &= 0 \\
# b:&& 2x + \phantom{2}y + z +0&= 0 \\
# c:&& x -2y -z -2 &= 0
# \end{aligned}
# $$

# %%
a = e1 + e2 + e3 - e0
b = 2*e1 + e2 + e3
c = e1 -2*e2-e3 +2*e0 

# %% [markdown]
# The solution is the intersection point of the 3 planes described by each equation above.

# %%
a^b^c

# %% [markdown]
# The trivector terms describe a point.  It's easier to read off if we normalize it and look at its orthogonal complement:

# %%
J3(a^b^c) 

# %% [markdown]
# The coefficient of $e_1$ is the $x$ value, etc.  So the solution of the original system of equations is
# $$
# \begin{aligned}
# x&=-1\\
# y&=-5\\
# z&=3
# \end{aligned}
# $$

# %%
# Let's define p as a generic point
p=J3(x*e1 + y*e2 + z*e3 + e0)
p

# %% [markdown]
# We can represent this as equations for x,y,z by setting terms of the expression below equal to zero: $p\vee(a\wedge b\wedge c)=0$.  PGA is a Regressive Product Null Space (RPNS) representation.

# %%
vee3(p,a^b^c)

# %% [markdown]
# Setting the line above equal to zero gives a redundant set of equations, but looking at the line three terms we can read off
# $$
# \begin{aligned}
# x&=-1\\
# y&=-5\\
# z&=7
# \end{aligned}
# $$
#
# ## Under-defined system
# Let's try something different.  Let's do an under-defined system.
# $$
# \begin{aligned}
# a:&& x + y + z -1 &= 0 \\
# b:&& 2x + y + z +0&= 0 \\
# %x -2y -z &= 2
# \end{aligned}
# $$
#
# Again, the solution is $a\wedge b$:

# %%
a^b

# %% [markdown]
# This is the Plücker representation of a line in 3D space.  As before, we can get this into a more familiar form by computing $p \vee (a\wedge b)=0$.

# %%
vee3(p,a^b)

# %% [markdown]
# ### Parametric form of the line
# (TODO generalize the method to higher dimensional solution spaces.)
#
# The solution to the system above is a line of points: $x=-1, y=-z+2$.
#
# To put this in parametric form, we need the direction of the line $d$ and any point on the line $P$:
# $$
# \begin{aligned}
# d &= a\wedge b \wedge e_0 \\
# P &= d^* \wedge (a\wedge b)
# \end{aligned}
# $$
#
# Then the line consists of the points $P + dt$ for a scalar value $t$.

# %%
d=a^b^e0 # direction of the line
d

# %%
P = J3(d)^(a^b) # a point on the line
P

# %%
(P+d*t)

# %% [markdown]
# Or, translating to coordinates and rescaling $t$:
# $$
# x = -1,\quad
# y = 1-t, \quad z=1+t
# $$
#
# ## Over-defined system
# Let's look at an over-defined system in 2D:
# $$
# \begin{aligned}
# x&= 0\\
# y&= 0\\
# x+y&= 1
# \end{aligned}
# $$
# This clearly has no solution, but can we find a "best fit"?
# One way to answer is to sum the 3 normalized intersection points.
#
# (Q: what happens if I use the _unnormalized_ points? Is this useful?)

# %%
sum([e1^e2,-e2^(e1+e2-e0),-(e1+e2-e0)^e1])/3

# %% [markdown]
# To generalize to larger systems, for an _N_-dimensional space with _m_ equations, form all $\binom{m}{N}$ of the _N_-wise wedge products, normalize, and sum them.
#
# This method also works when two or more of the (hyper-)planes are parallel.  The intersection point is then ideal, and that point can be discarded before summing.
#
# (Q: does this give the same result as the least-squares solution, using pseudo-inverse?  It seems like that is different?)
#
# ## Gaussian elimination
# One of the standard solution techniques in linear algebra is Gaussian elimination.  What does that look like in PGA?  It's almost trivial.
#
# Let's go back to our earlier problem:
# $$
# \begin{aligned}
# a:&&x + \phantom{2}y + z -1 &= 0 \\
# b:&&2x + \phantom{2}y + z +0 &= 0 \\
# c:&&x -2y -z -2 &= 0
# \end{aligned}
# $$
#
# The idea of Gaussian elimination is to add a scalar multiple of one row to another row to zero out coefficients in the lower left triangle of the matrix.
#
# For example, we could replace row $b$ with $b-2a$ and row $c$ with $c-a$, to eliminate the coefficients of $x$ in the bottom two rows:

# %%
b-2*a, c-a

# %% [markdown]
# $$
# \begin{aligned}
# a:&&x + \phantom{3}y +\phantom{2} z -1 &= 0 \\
# b-2a:&&0x  -\phantom{3}y - \phantom{2}z +2&= 0 \\
# c-\phantom{2}a:&&0x -3y -2z -1 &= 0
# \end{aligned}
# $$
#
# Then add -3 times the second row to the third row:

# %%
-3*(b-2*a)+(c-a)

# %% [markdown]
# $$
# \begin{aligned}
# a:&&x + y + \phantom{0}z -1 &= 0 \\
# b-2a:&&0x  -y - \phantom{0}z +2&= 0 \\
# 5a-3b+c:&&0x+ 0y +z -3 &= 0
# \end{aligned}
# $$
#
# But the properties of the exterior product make this unnecessary.
# $$
# a\wedge(b-2a)\wedge(5a-3b+c) = a\wedge b \wedge c,
# $$
# so we are back where we started.
#
# ### Gaussian elimination as a series of rotations
# We can also take a more geometric approach by thinking of Gaussian elimination as a series of rotations of the planes about their common lines of intersections.

# %%
test=normalize(1+normalize(a)*normalize(b))

# %%
X=e1; Y=e2; Z=e3

# %% [markdown]
# So, each step can take the form (for example)
# $
# b' = (e_1\cdot b)a
# $.  How do I gradually turn this on?

# %%
# Step 1 is this rotor
test=normalize(1+normalize(b)*normalize(X|b*a))
# I need the log of this.  It's a rotation, so should be of the form
# cos(a/2)+sin(a/2)B

# %%
angle=acos(test[0].obj)
B = test[2]/sin(angle)
cos(angle)+sin(angle)*B

# %%
#Step 2
c-a, X|c*a

# %%
#Step 3
Y|(X|c)

# %% [markdown]
# ### TODOs
# * add figures, highlight geometry
# * describe rotation picture of Gaussian elimination
# * projections, rejections, reflections
# * emphasize linear combinations of planes keep their intersection invariant
# * add norm and normalize functions
# * use macros to make the text automatically update if the equations are changed
