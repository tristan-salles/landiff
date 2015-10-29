# landiff
Multi-lithology non-linear diffusion model of hillslope dynamic

We consider here a sedimentation scenario of $n$ lithology. The following set of nonlinear PDEs (n+1 equations in total), derived from Riven\ae s, constitute the mathematical model:

$$\sum_{k=1}^n s_k = 1$$
 
$$
	\frac{\partial h}{\partial t} = \sum_{k=1}^n \frac{1}{C_k} \nabla \cdot \left( \kappa_k s_k \nabla h \right),
$$
and for each sediment $k$ in $\left[ 1,...,n-1 \right]$:
$$
	A \frac{\partial s_k}{\partial t} + s_k \frac{\partial h}{\partial t} = \frac{1}{C_k} \nabla \cdot \left( \kappa_k s_k \nabla h \right)
$$
