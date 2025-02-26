## Problem Description


The problem consists of simulating the heat diffusion equation:

```math
\frac{\partial \rho (x, y, t) }{\partial t} = D \cdot \nabla^2 \rho(x, y,t)
```

where $\rho$ is the temperature function (in Kelvin) and $D$ is the thermal diffusivity coefficient (in square meters per second).

We want to simulate the equation in a discretized square grid with side N.

By applying the Implicit Euler method we find:

```math
T^{n+1}_{x,y} = T^{n}_{x, y} + D \left( \frac{T^{n+1}_{x-1, y} + T^{n+1}_{x, y-1} - 4T^{n+1}_{x, y} + T^{n+1}_{x+1, y} + T^{n+1}_{x, y+1}}{\Delta h^2} \right) \Delta t
```

Once the terms are rearranged, we obtain:

```math
\frac{D \Delta t}{\Delta h^2} \left(4T^{n+1}_{x, y} - T^{n+1}_{x-1, y} - T^{n+1}_{x, y-1} - T^{n+1}_{x+1, y} - T^{n+1}_{x, y+1}\right) + T^{n+1}_{x, y}  = T^{n}_{x, y}
```

which holds in the case where there are no constraints at the internal grid points and if the points are not on the grid edges.

In the latter case, where the points are on the edges, we need to modify the equation accordingly.  
For example, in the case of the first point of the grid, where $(x, y) = (0,0)$, we have:

```math
\frac{D \Delta t}{\Delta h^2} \left(4T^{n+1}_{0, 0} - T^{n+1}_{1, 0} - T^{n+1}_{0, 1}\right) + T^{n+1}_{0, 0}  = T^{n}_{0, 0} + \frac{D \Delta t}{\Delta h^2} (C_{-1,0} + C_{0, -1})
```

where $C_{-1,0}$ and $C_{0,-1}$ represent the constant values (in Kelvin) of the cells just beyond the grid, and these values are fixed, they do not change with time.

Thus, if we only consider boundary conditions (immediately after the grid boundary),  
we can rewrite our problem in matrix form:

$$
A T^{n+1} = T^{n} + F
$$

where $A$ is the matrix containing the coefficients, $F$ contains any constant terms due to the boundary conditions, $T^n$ is the temperature vector at this time (we know it) and $T^{n+1}$ is the temperature vector at the next time step, the one we are interested in.

As an example, if the side length of our matrix is $N=4$, our initial grid will be:

```math
\left|
\begin{array}{c|cccc|c}
\hline
 & C_{-1, 0} & C_{-1, 1} & C_{-1, 2} & C_{-1, 3} & \\ \hline
C_{0, -1} & T_{00} & T_{01} & T_{02} & T_{03} & C_{0, 4}\\
C_{1, -1} & T_{10} & T_{11} & T_{12} & T_{13} & C_{1, 4}\\
C_{2, -1} & T_{20} & T_{21} & T_{22} & T_{23} & C_{2, 4}\\
C_{3, -1} & T_{30} & T_{31} & T_{32} & T_{33} & C_{3, 4}\\ \hline
& C_{4, 0} & C_{4, 1} & C_{4, 2} & C_{4, 3} & \\ \hline
\end{array}
\right|
```

This leads to the linear system $A T^{n+1} = T^{n} + F$:

```math
\left[
\begin{array}{cccc | cccc | cccc | cccc}
\alpha&\beta&&          &\beta&&&               &&&&                &&&& \\
\beta&\alpha&\beta&     &&\beta&&               &&&&                &&&& \\
&\beta&\alpha&\beta     &&&\beta&               &&&&                &&&& \\
&&\beta&\alpha          &&&&\beta               &&&&                &&&& \\ \hline
\beta&&&&               \alpha&\beta&&          &\beta&&&            &&&& \\
&\beta&&&               \beta&\alpha&\beta&     &&\beta&&            &&&& \\
&&\beta&&               &\beta&\alpha&\beta     &&&\beta&            &&&& \\
&&&\beta&               &&\beta&\alpha          &&&&\beta            &&&& \\ \hline
&&&&                    \beta&&&&               \alpha&\beta&&          &\beta&&& \\
&&&&                    &\beta&&&               \beta&\alpha&\beta&     &&\beta&& \\
&&&&                    &&\beta&&               &\beta&\alpha&\beta     &&&\beta& \\
&&&&                    &&&\beta&               &&\beta&\alpha          &&&&\beta \\ \hline
&&&&        &&&&                    \beta&&&&               \alpha&\beta&&     \\
&&&&        &&&&                    &\beta&&&               \beta&\alpha&\beta&\\
&&&&        &&&&                    &&\beta&&               &\beta&\alpha&\beta\\
&&&&        &&&&                    &&&\beta&               &&\beta&\alpha     \\
\end{array}
\right]
\cdot 
\left[
\begin{array}{c}
T"_{00} \\ T"_{01} \\ T"_{02} \\ T"_{03} \\
T"_{10} \\ T"_{11} \\ T"_{12} \\ T"_{13} \\
T"_{20} \\ T"_{21} \\ T"_{22} \\ T"_{23} \\
T"_{30} \\ T"_{31} \\ T"_{32} \\ T"_{33}
\end{array}
\right]
=
\left[
\begin{array}{c}
T_{00} \\ T_{01} \\ T_{02} \\ T_{03} \\
T_{10} \\ T_{11} \\ T_{12} \\ T_{13} \\
T_{20} \\ T_{21} \\ T_{22} \\ T_{23} \\
T_{30} \\ T_{31} \\ T_{32} \\ T_{33}
\end{array}
\right]
+
\left[
\begin{array}{c}
-\beta (C_{0,-1} + C_{-1,0}) \\ -\beta C_{-1, 1} \\ -\beta C_{-1, 2} \\ -\beta (C_{-1,3} + C_{0,4})  \\
-\beta C_{1,-1} \\ 0 \\ 0 \\ -\beta C_{1,4} \\
-\beta C_{2,-1} \\ 0 \\ 0 \\ -\beta C_{2,4} \\
-\beta (C_{3,-1} + C_{4,0}) \\ -\beta C_{4, 1} \\ -\beta C_{4, 2} \\ -\beta (C_{4,3} + C_{3,4}) 
\end{array}
\right]
```
where $\alpha = \left(1 + 4\frac{D \Delta t}{\Delta h^2} \right)$ while $\beta = -\frac{D \Delta t}{\Delta h^2}$.
Instead, $F$ is the vector that includes the boundary conditions.
If this is zero, it is equivalent to implicitly assuming that the temperature outside the grid is 0 Kelvin.
We used $`T_{xy}`$ as short for $`T^n_{xy}`$ and $` T"_{xy}`$ as short for $`T^{n+1}_{xy} `$.

If we want to modify these boundary conditions, we simply need to change the various $C$ values and adjust the vector $F$ accordingly.

Additionally, we note that our matrix is symmetric, which implies that its inverse is also symmetric and can therefore be represented using half the space. (Although computing the inverse still requires at least the full space).

We need to solve the system of equations to find the unknowns corresponding to the grid temperatures at the next time step. The problem is that the number of elements in the matrix scales as $N^4$, and consequently, if we want to represent our matrix in memory, we would need (assuming 8 bytes per value) approximately $762$ MB for $N=100$ or approximately $12$ GB for $N=200$.

Clearly, storing such a matrix in memory does not make sense, especially because it is zero almost everywhere except at specific points. (The proportion of nonzero values is given by: $\frac{5n-4}{n^3}$).

In general, solving the system does not require finding the inverse, and it is usually not advisable to do so. However, in our case, since the matrix associated with our system does not change over time, it is beneficial to compute it once, reducing the problem to a simple matrix-vector multiplication to determine the temperatures at the next time step:  

```math
T^{n+1} = A^{-1} (T^n + F)
```

Moreover, storing the matrix efficiently in memory is not difficult, as the relevant information is contained only in the five diagonals, each at most of length $N$. The issue is that the inverse matrix, in general, is not as sparse or structured in a predictable pattern.  

The implemented solution aims to optimize both memory usage and computational time for advancing the time step, ensuring that neither resource is excessively consumed.  

In case of a tridiagonal matrix, **Thomas algorithm** can calculate the inverse in linear time. We do not have a tridiagonal matrix, but if we consider the blocks, we have a block tridiagonal matrix:

```math
\begin{bmatrix}
D_0 & U_0 &  & & \\
L_1 & D_1 & U_1 & \\
 &  \ddots & \ddots & \ddots & \\ 
&  &  L_{N-2} & D_{N-2} & U_{N-2}\\
 &  &  & L_{N-1} & D_{N-1}\\ 
\end{bmatrix}
\cdot
\begin{bmatrix}
    T"_{0} \\ T"_{1} \\ \vdots \\ T"_{N-2} \\ T"_{N-1}
\end{bmatrix}
 =
\begin{bmatrix}
    T_{0} \\ T_{1} \\ \vdots \\ T_{N-2} \\ T_{N-1}
\end{bmatrix}
+
\begin{bmatrix}
    F_{0} \\ F_{1} \\ \vdots \\ F_{N-2} \\ F_{N-1}
\end{bmatrix}
```

where the blocks $D$ are tridiagonal matrices, while $U$ and $L$ are diagonal matrices:

```math
    D_i = 
    \begin{bmatrix}
     \alpha & \beta & &  \\
    \beta & \alpha & \beta &  \\
    & \beta & \alpha & \beta  \\
    & & \beta & \alpha  
    \end{bmatrix}
    \quad
    U_i = L_i = 
    \begin{bmatrix}
     \beta & & &  \\
    & \beta & &  \\
    & & \beta &  \\
    & & & \beta  
    \end{bmatrix}
    \quad
    T_i = 
    \begin{bmatrix}
        T_{i0} \\ T_{i1} \\ T_{i2} \\ T_{i3} 
    \end{bmatrix}
```
We can then use the Thomas Block Algorithm (put the citation) to solve a block tridiagonal matrix:

*TODO PUT IMAGE OF THE ALGO*

Where the operations need to be thinked as between matrices and vectors.

The advantage is that the problem is then reduced to computing the inverse of $N$ square matrices of size $N \times N$ once, creating auxiliary $N \times N$ matrices that allow for computing the next step using only algebraic operations between matrices.



To be continued
