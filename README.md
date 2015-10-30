# landiff - a simple parallel multi-lithology non-linear diffusion model of hillslope dynamic

## The problem

The simple creep diffusion equation for single lithology has the basic form:

![equation](http://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%20h%7D%7B%5Cpartial%20t%7D%20%3D%20%5Cnabla%20%5Ccdot%20%5Cleft%28%20%5Ckappa%20%5Cnabla%20h%20%5Cright%29)

where *h(x, y, t)* is the height above some arbitrary horizontal surface in the *x-y* plane, *t* is time and ![equation](http://latex.codecogs.com/png.latex?%5Cinline%20%5Ckappa)*(x,y,t)* is the transport or diffusion coefficient giving the effectiveness of the diffusion.

Using a parameter for the fraction of a given sediment, Rivenaes added a second equation to calculate the ratio of two sediments in a layer of deposited material. The modified equations include *s(x, y, t)* and *1-s(x,y,t)* as the fraction of the two sediments. 

In particular, we consider here a sedimentation scenario of *n* lithology. The following set of nonlinear PDEs (n+1 equations in total), derived from Rivenaes, constitute the mathematical model:

![equation](http://latex.codecogs.com/png.latex?%5Csum_%7Bk%3D1%7D%5En%20s_k%20%3D%201%2C)

![equation](http://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%20h%7D%7B%5Cpartial%20t%7D%20%3D%20%5Csum_%7Bk%3D1%7D%5En%20%5Cfrac%7B1%7D%7BC_k%7D%20%5Cnabla%20%5Ccdot%20%5Cleft%28%20%5Ckappa_k%20s_k%20%5Cnabla%20h%20%5Cright%29%2C)

and for each sediment *k* in *[ 1,...,n-1]*:

![equation](http://latex.codecogs.com/png.latex?A%20%5Cfrac%7B%5Cpartial%20s_k%7D%7B%5Cpartial%20t%7D%20&plus;%20s_k%20%5Cfrac%7B%5Cpartial%20h%7D%7B%5Cpartial%20t%7D%20%3D%20%5Cfrac%7B1%7D%7BC_k%7D%20%5Cnabla%20%5Ccdot%20%5Cleft%28%20%5Ckappa_k%20s_k%20%5Cnabla%20h%20%5Cright%29)

In the above model, ![equation](http://latex.codecogs.com/png.latex?%5Cinline%20%5Ckappa_k%28x%2Cy%29)  denote the diffusion coefficients for sediment *k*. In addition, ![equation](http://latex.codecogs.com/png.latex?%5Cinline%20C_k) is the compaction ratio of sediment-*k* type. Moreover, *A* is a constant representing the thickness of a prescribed top layer, in which sediments are transported.

The initial conditions are of the form ![equation](http://latex.codecogs.com/png.latex?%5Cinline%20h%28x%2C%20y%2C%200%29%20%3D%20h%5E0%28x%2C%20y%29) and ![equation](http://latex.codecogs.com/png.latex?%5Cinline%20s%28x%2C%20y%2C%200%29%20%3D%20s%5E0%28x%2C%20y%29). As boundary conditions, most of the boundary has the no-flow condition, *i.e.*, the homogeneous Neumann boundary condition:

![equation](http://latex.codecogs.com/png.latex?%5Cinline%20%5Cpartial%20h/%5Cpartial%20n%20%3D%20%5Cpartial%20s/%5Cpartial%20n%20%3D%200)


## Temporal discretisation

We note ![equation](http://latex.codecogs.com/png.latex?%5Cinline%20%5CDelta%20t) the time step size. Let superscript *l* be the time level index, such that ![equation](http://latex.codecogs.com/png.latex?%5Cinline%20h%5El) denotes ![equation](http://latex.codecogs.com/png.latex?%5Cinline%20h%28x%2Cy%2Cl%5CDelta%20t%29) and ![equation](http://latex.codecogs.com/gif.latex?%5Cinline%20s_k%5El) denotes ![equation](http://latex.codecogs.com/gif.latex?%5Cinline%20s_k%28x%2Cy%2Cl%5CDelta%20t%29). Then, the temporal derivatives are simply approximated as

![equation](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%20h%7D%7B%5Cpartial%20t%7D%20%5Csimeq%20%5Cfrac%7Bh%5E%7Bl&plus;1%7D-h%5El%7D%7B%5CDelta%20t%7D%20%5C%2C%20%2C%20%5C%3B%20%5Cfrac%7B%5Cpartial%20s_k%7D%7B%5Cpartial%20t%7D%20%5Csimeq%20%5Cfrac%7Bs_k%5E%7Bl&plus;1%7D-s_k%5El%7D%7B%5CDelta%20t%7D)

The remaining task of temporal discretisation is to choose time level *l* or *l+1*, or a combination of both, at which the right-hand-side terms of the PDEs system are to be evaluated. Different strategies will give rise to **fully-explicit**, **semi-implicit** and **fully-implicit** schemes. Here we compute at the **fully-explicit** scheme.

### Fully-explicit scheme

To avoid solving systems of nonlinear algebraic equations, the right-hand-side terms of the PDEs system can use the already computed *h* and ![equation](http://latex.codecogs.com/gif.latex?%5Cinline%20s_k) values. More specifically, the equations are transformed as follows, by a fully-explicit temporal discretisation:

![equation](http://latex.codecogs.com/gif.latex?%5Csum_%7Bk%3D1%7D%5En%20s_k%5E%7Bl&plus;1%7D%20%3D%201)
![equation](http://latex.codecogs.com/gif.latex?%5Cfrac%7Bh%5E%7Bl&plus;1%7D-h%5El%7D%7B%5CDelta%20t%7D%20%3D%20%5Csum_%7Bk%3D1%7D%5En%20%5Cfrac%7B1%7D%7BC_k%7D%20%5Cnabla%20%5Ccdot%20%5Cleft%28%20%5Ckappa_k%20s%5El_k%20%5Cnabla%20h%5El%20%5Cright%29)

and for each sediment *k* in *[1,...,n-1]*:

![equation](http://latex.codecogs.com/gif.latex?A%20%5Cfrac%7Bs_k%5E%7Bl&plus;1%7D-s_k%5El%7D%7B%5CDelta%20t%7D%20&plus;%20s%5E%7Bl&plus;1%7D_k%20%5Cfrac%7Bh%5E%7Bl&plus;1%7D-h%5El%7D%7B%5CDelta%20t%7D%20%3D%20%5Cfrac%7B1%7D%7BC_k%7D%20%5Cnabla%20%5Ccdot%20%5Cleft%28%20%5Ckappa_k%20s%5El_k%20%5Cnabla%20h%5E%7Bl&plus;1%7D%20%5Cright%29)

It should be noted that *h* is to be updated before s during each time step. This is why the newly computed ![equation](http://latex.codecogs.com/gif.latex?%5Cinline%20h%5E%7Bl&plus;1%7D) (from the second equation) is immediately used to compute ![equation](http://latex.codecogs.com/gif.latex?%5Cinline%20s_k%5E%7Bl&plus;1%7D). Another remark is that ![equation](http://latex.codecogs.com/gif.latex?%5Cinline%20s_k%5E%7Bl&plus;1%7D), instead of ![equation](http://latex.codecogs.com/gif.latex?%5Cinline%20s_k%5E%7Bl), is used in the ![equation](http://latex.codecogs.com/gif.latex?%5Cinline%20s_k%5Cfrac%7B%5Cpartial%20h%7D%7B%5Cpartial%20t%7D) term on the left-hand side of the last equation. Numerical experiments show that this simple trick improves the numerical stability of this fully-explicit scheme, in which both ![equation](http://latex.codecogs.com/gif.latex?%5Cinline%20h%5E%7Bl&plus;1%7D) and ![equation](http://latex.codecogs.com/gif.latex?%5Cinline%20s_k%5E%7Bl&plus;1%7D) are computed straightforwardly. The scheme has **first-order accuracy in time**.

### Spatial discretisation

We choose finite differences to carry out the spatial discretisation. This is mostly motivated by the numerical and programming simplicity. It can be mentioned that other spatial discretisation techniques, such as finite elements, can also use the same temporal discretisation discussed above.

#### Treatment of diffusion

It is standard to use centred difference for the two diffusion terms on the right-hand side of equation 2, for obtaining second-order accuracy in space. For example, centred difference applied to the ![equation](http://latex.codecogs.com/gif.latex?%5Cinline%20%5Cnabla%20%5Ccdot%20%5Cleft%28%20%5Ckappa_k%20s_k%20%5Cnabla%20h%20%5Cright%29) term gives the following discretised form:

![equation](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Ckappa_%7Bk%5C%3B%20i&plus;%5Cfrac%7B1%7D%7B2%7D%2Cj%7D%20s_%7Bk%5C%3B%20%7Bi&plus;%5Cfrac%7B1%7D%7B2%7D%2Cj%7D%7D%20%28h_%7Bi&plus;1%2Cj%7D-h_%7Bi%2Cj%7D%29%20-%20%5Ckappa_%7Bk%5C%3B%20i-%5Cfrac%7B1%7D%7B2%7D%2Cj%7D%20s_%7Bk%5C%3B%20%7Bi-%5Cfrac%7B1%7D%7B2%7D%2Cj%7D%7D%20%28h_%7Bi%2Cj%7D-h_%7Bi-1%2Cj%7D%29%20%7D%7B%5CDelta%20x%5E2%7D)

![equation](http://latex.codecogs.com/gif.latex?&plus;%20%5Cfrac%7B%5Ckappa_%7Bk%5C%3B%20i%2Cj&plus;%5Cfrac%7B1%7D%7B2%7D%7D%20s_%7Bk%5C%3B%20%7Bi%2Cj&plus;%5Cfrac%7B1%7D%7B2%7D%7D%7D%20%28h_%7Bi%2Cj&plus;1%7D-h_%7Bi%2Cj%7D%29%20-%20%5Ckappa_%7Bk%5C%3B%20i%2Cj-%5Cfrac%7B1%7D%7B2%7D%7D%20s_%7Bk%5C%3B%20%7Bi%2Cj-%5Cfrac%7B1%7D%7B2%7D%7D%7D%20%28h_%7Bi%2Cj%7D-h_%7Bi%2Cj-1%7D%29%20%7D%7B%5CDelta%20y%5E2%7D%20%5C%3B%20%2C)

where the subscripts *i,j* are the mesh point index for a 2D uniform grid with mesh spacing ![equation](http://latex.codecogs.com/gif.latex?%5CDelta%20x) and ![equation](http://latex.codecogs.com/gif.latex?%5CDelta%20y). In the above formula, the half-indexed terms are to be evaluated as, *e.g.*, ![equation](http://latex.codecogs.com/gif.latex?%5Ckappa_%7Bk%5C%3B%20i&plus;%5Cfrac%7B1%7D%7B2%7D%2Cj%7D%20s_%7Bk%5C%3B%20%7Bi&plus;%5Cfrac%7B1%7D%7B2%7D%2Cj%7D%7D%3D%28%5Ckappa_%7Bk%5C%3B%20i%2Cj%7Ds_%7Bk%5C%3B%20i%2Cj%7D%20&plus;%20%5Ckappa_%7Bk%5C%3B%20i&plus;1%2Cj%7Ds_%7Bk%5C%3B%20i&plus;1%2Cj%7D%20%29/2).

#### Treatment of convection

The equation 3 is a convection equation with respect to ![equation](http://latex.codecogs.com/gif.latex?s_k), because of the ![equation](http://latex.codecogs.com/gif.latex?%5Cnabla%20%5Ccdot%20%5Cleft%28%20%5Ckappa_k%20s_k%20%5Cnabla%20h%20%5Cright%29) term. For the sake of numerical stability, one-sided upwind finite difference is preferred over centred difference, despite its first-order spatial accuracy.

To this end, it is customary to move the convection term ![equation](http://latex.codecogs.com/gif.latex?%5Cnabla%20%5Ccdot%20%5Cleft%28%20%5Ckappa_k%20s_k%20%5Cnabla%20h%20%5Cright%29) to the left-hand side of equation when checking the flow direction. That is, ![equation](http://latex.codecogs.com/gif.latex?-%5Cnabla%20h) gives the convection velocity. The *x*-component, ![equation](http://latex.codecogs.com/gif.latex?-%5Cpartial%20h/%5Cpartial%20x), is approximated by ![equation](http://latex.codecogs.com/gif.latex?h_%7Bi-1%2Cj%7D-h_%7Bi&plus;1%2Cj%7D%29/2%5CDelta%20x), the sign of which determines how the *x*-component of the convection term is discretised by one-sided upwind difference. More specifically, the

![equation](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20x%7D%20%5Cleft%28%20%5Ckappa_k%20s_k%20%5Cfrac%7B%5Cpartial%20h%7D%7B%5Cpartial%20x%7D%20%5Cright%29) term is approximated by

![equation](http://latex.codecogs.com/gif.latex?%5Cleft%28%20%5Cfrac%7B%5Ckappa_%7Bk%5C%3B%20%7Bi%2Cj%7D%7D%20s_%7Bk%5C%3B%20%7Bi%2Cj%7D%7D%20-%20%5Ckappa_%7Bk%5C%3B%20%7Bi-1%2Cj%7D%7D%20s_%7Bk%5C%3B%20%7Bi-1%2Cj%7D%7D%20%7D%7B%20%5CDelta%20x%7D%20%5Cright%29%20%5Ctimes%20%5Cleft%28%20%5Cfrac%7Bh_%7Bi&plus;1%2Cj%7D-h_%7Bi-1%2Cj%7D%7D%7B2%5CDelta%20x%7D%20%5Cright%29) if we have ![equation](http://latex.codecogs.com/gif.latex?h_%7Bi-1%2Cj%7D%20%3E%20h_%7Bi&plus;1%2Cj%7D). Otherwise, the following approximation is used:

![equation](http://latex.codecogs.com/gif.latex?%5Cleft%28%20%5Cfrac%7B%5Ckappa_%7Bk%5C%3B%20%7Bi&plus;1%2Cj%7D%7D%20s_%7Bk%5C%3B%20%7Bi&plus;1%2Cj%7D%7D%20-%20%5Ckappa_%7Bk%5C%3B%20%7Bi%2Cj%7D%7D%20s_%7Bk%5C%3B%20%7Bi%2Cj%7D%7D%20%7D%7B%20%5CDelta%20x%7D%20%5Cright%29%20%5Ctimes%20%5Cleft%28%20%5Cfrac%7Bh_%7Bi&plus;1%2Cj%7D-h_%7Bi-1%2Cj%7D%7D%7B2%5CDelta%20x%7D%20%5Cright%29)

### Treatment of boundary conditions

Second-order accurate treatment of the homogeneous Neumann condition ![equation](http://latex.codecogs.com/gif.latex?%5Cpartial%20h/%5Cpartial%20n%20%3D%20%5Cpartial%20s/%5Cpartial%20n%20%3D%200) follows the standard approach by using one layer of ghost boundary points.

## Numerical scheme

The complete numerical discretisation of the equation 2 in the PDEs is as follows:

![equation](
http://latex.codecogs.com/gif.latex?%5Cbegin%7Baligned%7D%20%5Cfrac%7Bh%5E%7Bl&plus;1%7D_%7Bi%2Cj%7D-h%5El_%7Bi%2Cj%7D%7D%7B%5CDelta%20t%7D%20%3D%20%5Csum%5En_%7Bk%3D1%7D%20%5Cleft%5B%20%5Cfrac%7B1%7D%7BC_k%7D%20%5Cfrac%7B%20%5Ckappa_%7Bk%5C%3B%20i&plus;%5Cfrac%7B1%7D%7B2%7D%2Cj%7D%20s%5El_%7Bk%5C%3B%20%7Bi&plus;%5Cfrac%7B1%7D%7B2%7D%2Cj%7D%7D%20%28h%5El_%7Bi&plus;1%2Cj%7D-h%5El_%7Bi%2Cj%7D%29%20-%20%5Ckappa_%7Bk%5C%3B%20i-%5Cfrac%7B1%7D%7B2%7D%2Cj%7D%20s%5El_%7Bk%5C%3B%20%7Bi-%5Cfrac%7B1%7D%7B2%7D%2Cj%7D%7D%20%28h%5El_%7Bi%2Cj%7D-h%5El_%7Bi-1%2Cj%7D%29%20%7D%7B%5CDelta%20x%5E2%7D%20%5Cright.%20%5C%5C%20&plus;%20%5Cleft.%20%5Cfrac%7B1%7D%7BC_k%7D%20%5Cfrac%7B%20%5Ckappa_%7Bk%5C%3B%20i%2Cj&plus;%5Cfrac%7B1%7D%7B2%7D%7D%20s%5El_%7Bk%5C%3B%20%7Bi%2Cj&plus;%5Cfrac%7B1%7D%7B2%7D%7D%7D%20%28h%5El_%7Bi%2Cj&plus;1%7D-h%5El_%7Bi%2Cj%7D%29%20-%20%5Ckappa_%7Bk%5C%3B%20i%2Cj-%5Cfrac%7B1%7D%7B2%7D%7D%20s%5El_%7Bk%5C%3B%20%7Bi%2Cj-%5Cfrac%7B1%7D%7B2%7D%7D%7D%20%28h%5El_%7Bi%2Cj%7D-h%5El_%7Bi%2Cj-1%7D%29%20%7D%7B%5CDelta%20y%5E2%7D%20%5Cright%5D%20%2C%20%5Cend%7Baligned%7D
)

where as previously defined, the half-index subscripts in the above formula mean some form of averaging for a quantity in the middle of two spatial mesh points. We typically use an arithmetic mean as follows:

![equation](
http://latex.codecogs.com/gif.latex?%5Ckappa_%7Bk%5C%3B%20i&plus;%5Cfrac%7B1%7D%7B2%7D%2Cj%7D%20s%5El_%7Bk%5C%3B%20%7Bi&plus;%5Cfrac%7B1%7D%7B2%7D%2Cj%7D%7D%20%3D%20%5Cfrac%7B%20%5Ckappa_%7Bk%5C%3B%20i%2Cj%7D%20s%5El_%7Bk%5C%3B%20i%2Cj%7D%20&plus;%20%5Ckappa_%7Bk%5C%3B%20i&plus;1%2Cj%7Ds%5El_%7Bk%5C%3B%20i&plus;1%2Cj%7D%20%7D%7B2%7D
)

For lithology *k*, equation 3 the complete numerical discretisation reads:

![equation](http://latex.codecogs.com/gif.latex?%5Csmall%20%5Cbegin%7Baligned%7D%20A%20%5Cfrac%7Bs_%7Bk%5C%3B%20i%2Cj%7D%5E%7Bl&plus;1%7D-s_%7Bk%5C%3B%20i%2Cj%7D%5El%7D%7B%5CDelta%20t%7D%20&plus;%20s%5E%7Bl&plus;1%7D_%7Bk%5C%3B%20i%2Cj%7D%20%5Cfrac%7Bh_%7Bi%2Cj%7D%5E%7Bl&plus;1%7D-h_%7Bi%2Cj%7D%5El%7D%7B%5CDelta%20t%7D%20%3D%20%5Cfrac%7B1%7D%7B2C_k%5CDelta%20x%5E2%7D%20%5Cbegin%7Bcases%7D%20%5Cleft%28%20%5Ckappa_%7Bk%5C%3B%20%7Bi%2Cj%7D%7D%20s%5El_%7Bk%5C%3B%20%7Bi%2Cj%7D%7D%20-%20%5Ckappa_%7Bk%5C%3B%20%7Bi-1%2Cj%7D%7D%20s%5El_%7Bk%5C%3B%20%7Bi-1%2Cj%7D%7D%20%5Cright%29%20%5Cleft%28%20h%5E%7Bl&plus;1%7D_%7Bi&plus;1%2Cj%7D-h%5E%7Bl&plus;1%7D_%7Bi-1%2Cj%7D%20%5Cright%29%20%5C%3B%20%28c1%29%20%5C%5C%20%5Cleft%28%20%5Ckappa_%7Bk%5C%3B%20%7Bi&plus;1%2Cj%7D%7D%20s%5El_%7Bk%5C%3B%20%7Bi&plus;1%2Cj%7D%7D%20-%20%5Ckappa_%7Bk%5C%3B%20%7Bi%2Cj%7D%7D%20s%5El_%7Bk%5C%3B%20%7Bi%2Cj%7D%7D%20%5Cright%29%20%5Cleft%28%20h%5E%7Bl&plus;1%7D_%7Bi&plus;1%2Cj%7D-h%5E%7Bl&plus;1%7D_%7Bi-1%2Cj%7D%20%5Cright%29%20%5C%3B%20%28c2%29%20%5Cend%7Bcases%7D%20%5C%5C%20&plus;%20%5Cfrac%7B1%7D%7B2C_k%5CDelta%20y%5E2%7D%20%5Cbegin%7Bcases%7D%20%5Cleft%28%20%5Ckappa_%7Bk%5C%3B%20%7Bi%2Cj%7D%7D%20s%5El_%7Bk%5C%3B%20%7Bi%2Cj%7D%7D%20-%20%5Ckappa_%7Bk%5C%3B%20%7Bi%2Cj-1%7D%7D%20s%5El_%7Bk%5C%3B%20%7Bi%2Cj-1%7D%7D%20%5Cright%29%20%5Cleft%28%20h%5E%7Bl&plus;1%7D_%7Bi%2Cj&plus;1%7D-h%5E%7Bl&plus;1%7D_%7Bi%2Cj-1%7D%20%5Cright%29%20%5C%3B%20%28c3%29%20%5C%5C%20%5Cleft%28%20%5Ckappa_%7Bk%5C%3B%20%7Bi%2Cj&plus;1%7D%7D%20s%5El_%7Bk%5C%3B%20%7Bi%2Cj&plus;1%7D%7D%20-%20%5Ckappa_%7Bk%5C%3B%20%7Bi%2Cj%7D%7D%20s%5El_%7Bk%5C%3B%20%7Bi%2Cj%7D%7D%20%5Cright%29%20%5Cleft%28%20h%5E%7Bl&plus;1%7D_%7Bi%2Cj&plus;1%7D-h%5E%7Bl&plus;1%7D_%7Bi%2Cj-1%7D%20%5Cright%29%20%5C%3B%20%28c4%29%20%5Cend%7Bcases%7D%20%5Cend%7Baligned%7D)


with condition *(c1)* corresponding to ![equation](http://latex.codecogs.com/gif.latex?%5Csmall%20h%5E%7Bl&plus;1%7D_%7Bi-1%2Cj%7D%20%3E%20h%5E%7Bl&plus;1%7D_%7Bi&plus;1%2Cj%7D) and *(c2)* otherwise. Similarly condition *(c3)* is used in case ![equation](http://latex.codecogs.com/gif.latex?%5Csmall%20h%5E%7Bl&plus;1%7D_%7Bi%2Cj-1%7D%20%3E%20h%5E%7Bl&plus;1%7D_%7Bi%2Cj&plus;1%7D) and *(c4)* otherwise.

## Parallelisation

The parallelisation is simply done by splitting the regular grid between processors using either the rows or columns.

## Examples

Two examples are provided for testing purposes. You will need to manually edit the main files to change some specific variables such as the number of lithology or the coefficient of diffusion for marine and aerial environment.

The results are CSV files that could be open in Paraview... or anything else!


## References

- S.R. Clark, W. Wei, X. Cai, 2010. *Numerical analysis of a dual-sediment transport model applied to Lake Okeechobee, Florida*. in: Proceedings of the 9th International Symposium on Parallel and Distributed Computing, IEEE Computer Society Press, pp. 189-194.

- J.C. Rivenaes, 1992. *A computer simulation model for siliclastic basin stratigraphy*. Ph.D. thesis, University of Trondheim.

- J.C. Rivenaes, 1997. *Application of a dual-lithology, depth-dependent diffusion equation in stratigraphic simulation*. Basin Research, 4 (2), 133-146.

- W. Wei, S.R. Clark, H. Su, M. Wen, X. Cai, 2012. *Balancing efficiency and accuracy for sediment transport simulations*, [link](http://heim.ifi.uio.no/xingca/Wei-etal-2012-CG.pdf).

- H. Su, N. Wu, M. Wen, C, Zhang, X. Cai, 2013. *Performance of Sediment Transport Simulations on NVIDIA's Kepler Architecture*, in: International Conference on Computational Science, ICCS, Procedia Computer Science 18, 1275--1281.
