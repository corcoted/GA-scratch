# Units in PGA (Ted's version)
## Two dimensions
The basis vectors in 2D PGA are $e_0, e_1, e_2$. The question is what dimensional units to assign to these, if any.
Some requirements or at least desirable properties
* The elements of any $k$-vector should have the same dimensionality.
* Functions of multivectors should stay dimensionally sane (e.g. $\exp(B), B\in \bigwedge^2V$ should be unitless).

My proposal is that the basis element should be dimensionless, and interpreted as pure directions, even the ideal elements.
This requires that all of the coefficients of a $k$-vector have the same units -- Euclidean and ideal coefficients carry the same units.

So, a point might be
$p = (1\mathrm{cm}) e_{20} + (2\mathrm{cm}) e_{01} + (1\mathrm{cm}) e_{12}$.

What is the meaning of the units in the ideal coordinates?

## Camera model of projective geometry
(Nothing here requires the projective split, but it helps me to check things if I go to a concrete set of coordinates.)

To make sense of this, I think about the camera model of projective geometry.
For example, 2D projective geometry is like the image of the real 3D world onto a screen inside a pinhole camera.
For the point $p = xe_{20} + ye_{01} + f e_{12}$, the coefficient of the $e_{12}$ basis bivector is the distance from the pinhole to the screen ($f$ = focal length), and the other two coefficients are the $x,y$ coordinates as measured on the screen.

![image](Pinhole.pdf)

### Normalization
When we normalize a PGA point, we divide by the $e_{12}$ coefficient _including its dimensional units_.
This makes the normalized $x,y$ coefficients unitless as well.
How do we interpret $x,y$ then?  They are the "angles" measured at the pinhole between the camera axis and the point on the screen.  More specifically, in Euclidean PGA $x,y$ are the tangent of the actual angles in the 3D embedding space.

### Lines
This idea works for lines as well.  A line in 2D PGA is something like
$L = a e_1 + b e_2 + c e_0$.
Normalized this becomes
$$
\hat{L} = \frac{a}{\sqrt{a^2 + b^2}}\,e_1 +  \frac{b}{\sqrt{a^2 + b^2}}\,e_2 +  \frac{c}{\sqrt{a^2 + b^2}}\,e_0.
$$
Here the normalized coefficients of $e_1$ and $e_2$ are the direction cosines between the line and the respective axes, so clearly these are dimensionless.
We usually think of the $e_0$ coefficient as the distance between the line and the origin, but here it appears dimesionless.  So, what is it?  Again, it is an "angle" with respect to the origin of the embedding space (pinhole of the camera).  We can again identify this with the tangent of the actual angle in the embedding space.
The distance $1/\sqrt{a^2 + b^2}$ is again the distance from the pinhole to the screen.

$a,b$, and I propose $c$ above should have units of inverse length.

![image](line.pdf)

### Duality
(I think Leo would argue that points don't have any dimensionality, nor planes.  Only differences between points have a length.  In other words, the length is induced by the projective split.)

This one has me stumped.  Coincident points and lines can be expressed
$$
L \wedge p = (ax + by + cf)I = 0
$$
or
$$ L \cdot p^\star = ax + by + cf = 0. $$
The zeros make it difficult to deduce dimensionality, but we at least see that the pairs $ax$, $by$, and $cf$ must have the same units.

Looking at the dual of an object with itself, we have
$$ A \wedge A^\star = I = e_0e_1e_2\cdots$$
If the basis elements are dimensionless, then $A$ and $A^\star$ must have reciprocal dimensions.

This somewhat makes sense.  For example, points have dimensions of length and hyperplanes have dimensions of 1/length (as is usually the case in physics).

This relationship does, for example, help to distinguish between join lines and meet lines in 3D.
Join lines (join of two points) have dimensions of length² and meet lines (meet of two planes) have dimensions of length¯².  The moral of the story is that we cannot deduce the dimensionality from the grade alone, even for purely geometric objects.

Continuing: a join plane in 3D has dimensions of length³, and is dual to a meet point with dimensions of length¯³.

So, is there any meaningful way to relate a meet point and a elementary point at the same location?

Inverses are also problematic.  A point in 3D has grade 3 and dimensions length¹.  Its inverse will have grade 3 and dimensions length¯¹. This seems strange, but I guess it is not problematic.

### Other geometries
The idea of the norm of the element indicating the distance from the pinhole to the screen also makes sense in the case of hyperbolic and elliptic geometries.  The only thing different in those cases is the shape of the screen.  The elliptic case is not too hard to picture because it can be mapped onto spherical geometry.  If the norm of the vector/bivector represents the radius of the sphere, then the other normalized coefficients become angles in the literal sense.

## Weighted points/lines
So, what to make of weighted $k$-vectors in this picture?  Because the normalized $k$-vectors are dimensionless, any physical variable has to be weighted at least by some reference unit.  For example, to go from a normalized position to a "real" position, we need to multiply by the (scalar) unit of measure.  So we have two competing normalizations: (1) absolute normalization takes everything back to dimensionless "angles" and dimensionally normalized where we multiply back through by a scalar unit-ful constant.  This back-and-forth is a little cumbersome, but that's my best picture so far.

Mass points seem to be ok, dimensionally:
$$ M = \sum_i m_i X_i $$

Velocity looks a little funny.  This is a bivector (in 2D), so we should have something like
$$
V \stackrel{?}{=} \dot{x} e_{20} + \dot{y} e_{01} + \omega e_{12},
$$
where $\omega$ is angular velocity rotating about the origin (in the projective plane), and the over dots are time derivatives.

I think the camera model saves us again.  If the object is rotation with angular velocity $\omega$ in the projective plane, then its velocity in the embedding space is $\omega f$, where $f$ is the distance from the origin of the embedding space to the projective plane.

So, the corrected velocity bivector is
$$
V = \dot{x} e_{20} + \dot{y} e_{01} + f\omega e_{12}.
$$

With this, the momentum becomes (for a point mass)
$$
L = M\vee V = m\bigl[ (fy\omega-f\dot{y}) e_1 + (f\dot{x}-f\omega x) e_2 + (x\dot{y}-y\dot{x}) e_0 \bigr]
$$
I think this is correct.  The dimensions of each term are $[ML^2T^{-1}]$ like we expect.

# Response to Leo's draft
(Response to Leo's "Physical Geometry" chapter draft, dated 2021.11.18)

In Leo's draft, the ideal plane is represented by $e$ which is $-e_0$ in my notation above.
(Why the choice of minus sign?)
Leo handles the units question by assigning a dimension of 1/length to $e$.
With this choice, points become dimensionless (unit-wise as well as geometrically).
Extending this, all ideal elements have an extra factor of 1/length compared to their Euclidean counterparts.
I think (but I need to check) that this will get confusing when we start looking at dual elements.

My objection is that I do not think the basis elements should have dimensional units attached.
This ties you to a particular measurement system.
My way of doing things does not apply units until a coordinate system is applied.

## Coordinate-free objects
I have no objection to anything Leo says here.
This is independent of the issue of assigning dimensional units.
Two sanity checks: First, as Leo says, no results should depend on the choice of origin.
Second, the coordinate-free expressions must be dimensionally consistent with what you get after applying a coordinate system.

## What is projective space?
I think one of the underlying assumptions is how we describe the projective space and the objects living within it, in contrast to "real" (Euclidean) objects.
For example, a real point is a dimensionless object.  We can only define its position in relation to another point (e.g. the origin).
But, what is a point when we map it onto the projective space?  Is it still a dimensionless object?
Geometrically, the point is a line through some (arbitrary) origin of the embedding space.
One question is whether this last statement is even necessary.
Is it necessary to have an embedding space, or is the embedding space just a convenient way to apply coordinates?
In pure projective geometry the only real quantitative measurement is the cross ratio between 4 co-linear points or coincident lines.
Projective space is not a metric space.  So, we're really dealing with something different in PGA.  We do have a metric, albeit a degenerate one.
I suppose nothing should depend on the choice of metric, so long as it has the correct signature (e.g. a nonorthogonal basis or curvilinear basis should give the same results in the end).

Tangent: how curvilinear coordinate systems work.

## Building from planes


## Duals, reciprocals, and normalization