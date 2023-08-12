# ---
# title: 'R(2,2) playground'
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-autoscroll,-collapsed,-scrolled,-trusted,-ExecuteTime
#     notebook_metadata_filter: kernelspec,jupytext,-jupytext.text_representation.jupytext_version
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python [conda env:ga]
#     language: python
#     name: conda-env-ga-py
# ---

# %% [markdown]
# Ted Corcovilos, 2021-01-08
#
#
# Playing around with the 2d "mother algebra" $R(2,2)$, as described in C. Doran, et al., "Lie Groups as Spin Groups," *Journal of Mathematical Physics 34*, 3642 (1993). doi:[10.1063/1.530050](http://doi.org/10.1063/1.530050)
#
# I'll name the basis vectors $p_1, p_2, m_1, m_2$ with the diagonal metric $[1,1,-1,-1]$.
#
# Position vectors in this basis are represented by null vectors of the form $x (p_1+m_1) + y (p_2+m_2)$.

# %%
from sympy import *
from galgebra.ga import Ga
from galgebra.printer import latex
from IPython.display import Math

init_printing(latex_printer=latex, use_latex='mathjax')

# %%
xy = (x,y) = symbols("x y", real=True)

# %%
xyxy = (xp, yp, xm, ym) = symbols("xp yp xm ym", real=True)

# %%
R22 = Ga('p1 p2 m1 m2', g=[1,1,-1,-1], coords=xyxy)

# %%
p1, p2, m1, m2 = R22.mv() # break out basis vectors

# %%
# a real position vector has the form...
r = x*(p1+m1)+y*(p2+m2)

# %%
r

# %%
a, b, c, d = symbols("a b c d", real=True)
f, g, h, j = symbols("f g h j", real=True)

# %%
# big pseudoscalar
I=p1^p2^m1^m2
# special bivectors (ref. Doran)
K = (p1^m1)+(p2^m2)
E = (p1^p2)+(m1^m2)
F = (p1^m2)+(p2^m1)

# %% [markdown]
# From the paper, only bivectors that commute with $K$ will preserve null vectors. $E$ corresponds to the pseudo-scalar in the 2d position space. $F$ is the leftover bit needed to complete the bivector space.  Geometrically, $K$ describes scaling, $E$ rotations, and $F$ shears.

# %%
B1 = a*(p1^m1)+b*(p2^m2)+c*(p1^m2)+d*(p2^m1)+f*(p1^p2)+g*(m1^m2)

# %%
# check if commutator is zero
B1 >> K

# %%
# So, this will be zero iff c=d and f=-g
# Redefine, change up the names a bit:
B1 = (a+b)/2*(p1^m1)+(a-b)/2*(p2^m2)+c/2*(p1^m2)+c/2*(p2^m1)+d/2*(p1^p2)-d/2*(m1^m2)

# %%
# check commutator again
B1 >> K

# %%
# also need to normalize B1
# the norm squared is
Bnorm2=(B1*(B1.rev())).scalar()

# %%
A, B, C, D = symbols("A B C D", real=True)

# %% [markdown]
# Let's break `B1` down term by term to see how it transforms position.

# %%
# Define a small number as a placeholder to go from members of the Lie algebra to the Lie group
ϵ = symbols("ϵ", real=True)

# %% [markdown]
# Look at each bivector in the Lie algebra and see how the corresponding Lie group member transforms a position vector.

# %%
# for the "a" term
((-ϵ*p1^m1/2).exp())*r*((ϵ*p1^m1/2).exp())

# %% [markdown]
# So, the $p_1 \wedge m_1$ term looks like a scaling of $x$.

# %%
# for the "b" term
((-ϵ*p2^m2/2).exp())*r*((ϵ*p2^m2/2).exp())

# %% [markdown]
# Scaling of $y$, as expected.

# %%
# confirm that the a and b pieces commute:
(p1^m1) >> (p2^m2)

# %%
# commuting Lie algebra elements => we can apply the exponentials individually
# overall scaling
(-ϵ*p1^m1/2).exp()*((-ϵ*p2^m2/2).exp())*r*((ϵ*p2^m2/2).exp())*(ϵ*p1^m1/2).exp()

# %%
# inverse scaling (is there a better name?)
(-ϵ*p1^m1/2).exp()*((ϵ*p2^m2/2).exp())*r*((-ϵ*p2^m2/2).exp())*(ϵ*p1^m1/2).exp()

# %%
# check the pieces of the c term for commuting
(p1^m2) >> (p2^m1)

# %%
# the c term 
(-ϵ*p1^m2/2).exp()*((-ϵ*p2^m1/2).exp())*r*((ϵ*p2^m1/2).exp())*(ϵ*p1^m2/2).exp()

# %% [markdown]
# Scissor shear? (Not really a boost because this isn't Minkowski space...)

# %%
#check the d term for commuting:
(p1^p2) >> (m1^m2)

# %%
# the d term
(-ϵ*p1^p2/2).exp()*(ϵ*m1^m2/2).exp()*r*(-ϵ*m1^m2/2).exp()*(ϵ*p1^p2/2).exp()

# %% [markdown]
# Rotation.
#
# So, the matrix version of these terms looks something like
# $$
# \begin{array}{cc}
# a \rightarrow \begin{pmatrix} e^A & 0 \\ 0 & e^A \end{pmatrix}
# &
# b \rightarrow \begin{pmatrix} e^B & 0 \\ 0 & e^{-B} \end{pmatrix}
# \\
# c \rightarrow \begin{pmatrix} \cosh C & \sinh C \\ \sinh C & \cosh C \end{pmatrix}
# &
# d \rightarrow \begin{pmatrix} \cos D & -\sin D \\ \sin D & \cos D \end{pmatrix}
# \end{array}
# $$
#
# Note that these do not mutually commute, so it's hard to decouple a generic matrix into these pieces.
#
# Also, these are all positive-definite matrices, so we're not covering the full GL group.  Need one more piece:
# $$
# m \rightarrow \begin{pmatrix} 1 & 0 \\ 0 & \pm 1 \end{pmatrix}
# $$
# to cover reflections.
#
# Stepping back, we can identify the 4 Lie algebra generators as the identity and the (real) Pauli matrices, demonstrating the isomorphism between the bivectors and spinors.

# %%
# for example σ_x:
simplify(exp(C*Matrix([[0,1],[1,0]])))

# %%
# σ_y
simplify(exp(D*Matrix([[0,-1],[1,0]])))

# %%
# σ_z
simplify(exp(B*Matrix([[1,0],[0,-1]])))

# %%
# I
simplify(exp(A*Matrix([[1,0],[0,1]])))

# %% [markdown]
# We can combine the last three terms and simplify using the axis-angle formula for SU(2).  (Need to use complex values?)
#
# This still leaves the open problem of decomposing a generic 2x2 matrix into a product of the operators above.   The linear algebra solution is straight-forward but tedious.  Does GA give us any shortcuts?
