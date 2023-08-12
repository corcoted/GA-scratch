# ---
# title: Statics
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
# Ted Corcovilos - 20210609
#
# Statics problems in PGA, based on the ideas in Gunn's thesis and preprint.

# %%
from clifford.pga import *

# %%
from math import *

# %%
layout

# %%
I = e0123 # shortcut for the pseudoscalar

# %% [markdown]
# ## Forces
# In PGA, forces and torques are represented by lines through the point the force is applied to.  If the force has magnitude $F$, direction $\mathbf{u}$ (3-vector) and is applied at point $\mathbf{r}$ (3-vector), then the force bivector $\bar{F}$ (in 3d) is
#
# $$
# \bar{F} = F \mathbf{r}\vee\mathbf{u}
# $$
#
# For example, let's calculate the gravity force on a particle with mass $m$ located at coordinates $(b,0,0)$ under gravitational acceleration in the $-\hat{z}$.

# %%
# define our constant values
b = 1.0; g = 10.0; m=1.0

# %%
r = (e0+b*e1).dual()

# %%
u = (-1*e3).dual()

# %%
F = m*g*(r & u) # & is the vee product

# %%
F

# %% [markdown]
# The acceleration spear is represented by
# $$ \bar{a} = \mathbf{r} \vee \ddot{\mathbf{r}} $$

# %%
r & (g*u)

# %%
(~r)*r

# %% [markdown]
# ## Example 1
# Mass supported by two strings.
#
# A 1-kg mass is supported by two stings each of length 5 cm.  The ends of the strings are tied to 2 posts separated by 8 cm and at the same height.  What is the tension force on each string?
#
# TODO insert picture
#
# We know from symmetry that $T_1 = T_2$, but let's ignore that for now.
#
# We'll place the mass at the origin.  Its position will just be the bivector (wlrking with 2DPGA)
# $$
# \bar{M} = e_{12}.
# $$
#
# The two posts will be located at positions
# $$
# \bar{P}_1 = (e_{12} - 4 e_{20} + 3 e_{01})\,\mathrm{cm}
# \qquad
# \bar{P}_2 = (e_{12} + 4 e_{20} + 3 e_{01})\,\mathrm{cm}
# $$
#
# The two tension forces will be along the lines joining the mass with the points.
# $$
# \begin{aligned}
# \bar{T}_1 &= T_1 \widehat{\bar{M} \vee \bar{P_1}} \\
#           &= T_1 \widehat{e_{12} \vee (e_{12} - 4e_{20} + 3e_{e01})} \\
#           &= T_1 \widehat{J[e_0 \wedge (e_0 - 4e_1 + 3e_2)]} \\
#           &= T_1 \widehat{J[-4e_{01}+3e_{02}]} \\
#           &= T_1 \widehat{-4e_2 -3e_1} \\
# \bar{T}_1 &= T_1 \left(-\frac35 e_1 - \frac45 e_2 \right)\\[1.5ex]
# \bar{T}_2 &= T_2 \left(-\frac35 e_1 + \frac45 e_2 \right)
# \end{aligned}
# $$
#
# Above, $\hat{x}$ means normalized $x$:
# $$ \hat{x} = \frac{x}{\lVert x \rVert}, $$
# and $J$ represents the dual operation.
#
# #### The Code
# Right now the code is set up with 3DPGA, so all of the expressions will look slightly different.

# %%
# magnitude of mass
m = 1
# gravitational acceleration (approx.)
g = 10

# %%
# position of the mass
# I'll define the positions "dually"
M = (e0).dual()
print("M=",M)

# %%
# positions of the posts (units of cm)
P1 = (e0 - 4*e1 + 3*e2).dual()
P2 = (e0 + 4*e1 + 3*e2).dual()
print("P1=",P1, "\nP2=",P2)

# %%
# line of tension T1 (unitless)
l1 = (M & P1).normal()
print("l1=",l1)

# %%
# line of tension T2 (units of N)
l2 = (M & P2).normal()
print("l2=",l2)

# %% [markdown]
# We'll express the gravity force using the recipe
# $$
# \textrm{Force} = \textrm{magnitude}\; (\textrm{position} \vee \textrm{direction}).
# $$

# %%
# gravity force (units of N)
Fg = m*g*(M & ((-e2).dual()))
print("Fg=", Fg)

# %% [markdown]
# In static equilibrium, the sum of all of the force multivectors must be zero.  We can do this for each component separately:
# First, along the $x$-axis ($e_{23}$)
# $$ e_{23}: 0 = T_1 (-0.8) + T_2 (+0.8) \quad\rightarrow\quad T_1 = T_2.$$
#
# Then along the $y$-axis ($e_{13}$):
# $$
# e_{13}: 0 = T_1(-0.6) + T_2 (-0.6) + 10.0\,\mathrm{N} \qquad\rightarrow\qquad T_1 = T_2 = 8.33\,\mathrm{N}.
# $$
#
# # Example 2
#
# A person lifts up and holds one end of a long plank of mass $M$ and length $d$.  The other end is fixed on the ground.  What is the force exerted by the person?
#
# Call the gravity force $F_g$, the person force $F$, and the normal force on the stationary end $F_N$.
# We'll pick the coordinate origin to be the stationary end of the plank.
#
# $$
# \begin{aligned}
# \bar{F}_g &= Mg \left( e_{12} + \frac{d}{2}\,e_{20} \right) \vee (-e_{01}) \\
#           &= Mg \left(-\frac{d}{2} e_0 + e_1 \right)\\[1.5ex]
# \bar{F}_N &= F_N e_{12} \vee e_{01} = -F_N e_1 \\[1.5ex]
# \bar{F}   &= F (e_{12} + d\,e_{20}) \vee e_{01} \\
#           &= -F e_1 + Fd e_0
# \end{aligned}
# $$
#
# As before, the sum of the forces must be zero.  Looking at the $e_0$ terms, we get
# $$
# e_0:\qquad -\frac12 Mgd + Fd = 0 \quad \rightarrow \quad
# F = \frac12 Mg
# $$
#
# # Example 3
# (insert figure)
#
# A triangular truss consists of 4 points, located at points (in meters)
# $$
# \begin{array}{cc}
# A = (0,0) & B = (1,0) \\
# C = (1/2, \sqrt{3}/2) & D = (3/2, \sqrt{3}/2)
# \end{array}
# $$
#
# Five rigid links connect the points: $AB$, $AC$, $BC$, $BD$, $CD$.  Joints $A$ and $B$ are fastened to the ground.
#
# A mass with weight $W = 100\,\mathrm{N}$ is attached to point $D$.
#
# What is the tensile force in each link of the truss?

# %%
# the points
A = (e0).dual()
B = (e0 + e1).dual()
C = (e0 + 0.5*e1 + sqrt(3)/2.*e2).dual()
D = (e0 + 1.5*e1 + sqrt(3)/2.*e2).dual()

# %%
# the lines of the links, normalized
AB = (A & B).normal()
AC = (A & C).normal()
BC = (B & C).normal()
BD = (B & D).normal()
CD = (C & D).normal()

# %%
# weight force (100 N at point D in the -y direction)
W = 100.*(D & (-e2.dual())).normal()

# %%
W
