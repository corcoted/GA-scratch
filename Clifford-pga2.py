# ---
# title: Testing 2d pga in `clifford`
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

# %%
from clifford import *

# %%
import numpy as np

# %%
# define 2D PGA
layout, blades = Cl(2,0,1, firstIdx=0)
locals().update(blades)

# %%
# define some aliases for the duals of the basis vectors
E1=e1.dual()
E2=e2.dual()
E0=e0.dual()

# %%
blades


# %%
def rand_line(num=1,length=1.0):
    '''
    generates 'num' random (normalized) lines with maximum moment of 'length'
    '''
    # a line looks like xe1 + ye2 + ce0 and is normalized if x**2 + y**2 = 1
    # So, we'll use a polar form for x and y
    angle = 2.0*np.pi* np.random.random(num) # random angle between 0 and 2pi radians
    x = np.cos(angle)
    y = np.sin(angle)
    c = np.random.random(num)*length
    if num == 1:
        line = MultiVector(layout,[0,c[0],x[0],y[0],0,0,0,0])
    else:
        line = [MultiVector(layout,[0,c[i],x[i],y[i],0,0,0,0]) for i in range(num)]
    return(line)


# %%
def rand_point(num=1,length=1.0):
    '''
    generates 'num' random (normalized) lines with maximum moment of 'length'
    '''
    # a normalized point looks like -xe02 + ye01 + e12 
    # So, we'll use a polar form for x and y
    angle = 2.0*np.pi* np.random.random(num) # random angle between 0 and 2pi radians
    c = np.random.random(num)*length
    x = np.cos(angle)*c
    y = np.sin(angle)*c
    if num == 1:
        line = MultiVector(layout,[0,0,0,0,y[0],-x[0],1.0,0])
    else:
        line = [MultiVector(layout,[0,0,0,0,y[i],-x[i],1.0,0]) for i in range(num)]
    return(line)


# %%
# TODO define euclidean norm and ideal norm
# TODO write functions to denerate ideal points and the ideal line

# %%
# some examples

# %%
# generate 3 points
A, B, C = rand_point(3,5.0)

# %%
# calculate (unnormalized) lines between the points
a = B.vee(C)
b = C.vee(A)
c = A.vee(B)
[a,b,c]

# %%
# ERROR the meet and join methods don't work

# %%
# Verify that the distance between the points is the magnitude of the unnormalized line
[((B-C).as_array()**2).sum() == a.mag2(),
((C-A).as_array()**2).sum() == b.mag2(),
((A-B).as_array()**2).sum() == c.mag2()]

# %%
# normalize the lines
ahat = a/np.sqrt(a.mag2())
bhat = b/np.sqrt(b.mag2())
chat = c/np.sqrt(c.mag2())
[ahat,bhat,chat]

# %%
# TODO load pyganja and draw this

# %%
# check that the lines intersect back at the corresponding points
[
(a^b)/np.sqrt((a^b).mag2()) == C ,
(b^c)/np.sqrt((b^c).mag2()) == A ,
(c^a)/np.sqrt((c^a).mag2()) == B
]

# %%
# Calculate the area of the triangle
0.5*a^b^c

# %% [markdown]
# ## Drawing
# Try to draw these with `pyganja`

# %%
import pyganja

# %%
scene = pyganja.GanjaScene()
for p in [a,b,c]:
    scene.add_object(p)
for L in [A,B,C]:
    scene.add_object(L)

# %%
pyganja.draw(scene, sig=layout.sig)

# %% [markdown]
# ## Ray tracing: thin lens
# The first problem we'll tackle is a thin lens with focal length _f_.  We'll trace some random rays through the lens by two methods:
# 1. ABCD matrix
# 2. A pair of reflections
#
# We'll work in the paraxial approximation for now.
#
# A ray with height _h_ and slope _m_ hitting a thin lens of focal length _f_ located at the origin will leave the lens with height _h'_ and slope _m'_ according to
# $$
# \begin{pmatrix} h' \\ m' \end{pmatrix}=
# \begin{pmatrix} 1 & 0 \\ -1/f & 1 \end{pmatrix}
# \begin{pmatrix} h \\ m \end{pmatrix}
# = \begin{pmatrix} h \\ m - (h/f) \end{pmatrix}
# $$

# %%
f = 2. # focal length
A = 1.; B = 0.; C = -1./f; D=1. # thin lens

# %%
# let's randomly generate some heights and slopes.  I'll pick this to be in the range -1..1
num = 5 # number of rays
h = np.random.random(num)*2.0 - 1.0
m = np.random.random(num)*2.0 - 1.0

# %%
# calculate the output rays' heights and slopes
hp = h
mp = m - h/f

# %%
# generate PGA lines corresponding to the rays
r = [(m[i]*e1 - 1.*e2 + h[i]*e0) for i in range(num)]

# %% tags=[]
r

# %%
#manually calculate the output rays
rp = [(mp[i]*e1 - 1.*e2 + hp[i]*e0) for i in range(num)]

# %% [markdown]
# We need to translate the "rays" of geometric optics into line objects in PGA.
# The line $a e_1 + b e_2 + c e_0$ corresponds to the line with equation $ax+by+c=0$.
# So, a ray with height _h_ (at the y axis) and slope _m_ obeys the equation $y=mx+h \rightarrow mx - y + h = 0$.
# This then becomes the line object $\ell=me_1 -1e_2 + he_0$.
#
# Our ray transfer matrix needs to be converted also:
# $$
# \begin{pmatrix} A & B \\ C & D \end{pmatrix}\begin{pmatrix} h \\ m \end{pmatrix}
# \rightarrow
# \begin{pmatrix} A & B & 0 \\ C & D & 0 \\ 0 & 0 & 1\end{pmatrix}
# \begin{pmatrix} h \\ m \\ -1 \end{pmatrix}
# $$

# %%
# create matrix and then generate the outermorphism
ABCD = np.array([[A,B,0.],[C,D,0.],[0.,0.,1.]])

# %%
from clifford import transformations

# %%
ABCD_OM = transformations.OutermorphismMatrix(ABCD,layout)

# %%
# check that it works
[ABCD_OM(r[i])==rp[i] for i in range(num)]

# %%
# Now apply it to a point and see what image comes out.
test_obj = (r[0]^r[1])

# %%
test_obj

# %%
test_obj=test_obj/test_obj.value[6] # normalize by dividing by e12 coefficient

# %%
# output (image point)
test_img=ABCD_OM(test_obj).normal()

# %%
test_img

# %% [markdown]
# We'll compare this to the result from Gauss's Thin Lens Equation:
# $$
# \frac{1}{s}+\frac{1}{s'} = \frac{1}{f} \rightarrow s' = \frac{sf}{s-f}
# $$

# %% tags=[]
# manually calculate the image distance and height, with sign flip for s'
s = test_obj.value[4] # pull out x coordinate
sp = -s*f/(s-f)
# compare the two values
[test_img.value[4],sp]

# %%
# now check the heights
ht = test_obj.value[5] # pull out y coordinate
htp = sp/s*ht
# compare the values
[test_img.value[5],htp]

# %% [markdown]
# ## TODO
# Plot this stuff.  Currently `pyganja` only supports 3D PGA, 2D CGA, and 3D PGA.  So, I'll need to write conversion functions for 2D PGA to 2D CGA objects to allow plotting.
