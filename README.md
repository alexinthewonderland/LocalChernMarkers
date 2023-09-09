# LocalChernMarkers
This page is currently is still on progress sorry! It is mainly a repository of my work done on calculating the Local Chern Marker (LCM) which is done under the supervision of Prof. Ady Stern from the Weizmann Institute of Science. 

In this repo, you will see a file (still also on progress) about what I learned and summary of this project. This note will be a continuously updated note, until the project is finished. 

# Important Files things to note
* TM_calculation.m
* LCM
* n_orbital_calculation



# Introduction 
## The Hamiltonian Model
We implement the following Hamiltonian model in the momentum space 

<p align="center">
$\ H(\textbf{k}) = d_z(\textbf{k}) \sigma_z \otimes I_2 + d_x(\textbf{k}) \sigma_x \otimes I_2 + d_y(\textbf{k}) \sigma_y \otimes I_2$
</p>

where $\ \sigma_i$ are pauli matrices and the coefficients in front of them are defined as

<p align="center">
  $\ d_z = m+t_1 \sin^2{(k_x/2)}+t_1 \sin^2{(k_y/2)}$
</p>
  
<p align="center">
  $\ d_x = t_2 \sin{(k_x)}$
</p>
  
<p align="center">
  $\ d_y = t_2 \sin{(k_y)}$
</p>

Hence, in matrix form, the Hamiltonian in momentum space would be

```math
H =
\sum_{k}
\bf{c_k}
\begin{bmatrix}d_z & 0 & d_x-id_y & 0 \\
0 & d_z & 0 & d_x-id_y \\
d_x+id_y & 0 & -d_z & 0 \\
0 & d_x+id_y & 0 & -d_z
\end{bmatrix}
\bf{c_k^\dagger}
```

in which $\ a$, $\ b$, $\ c$, and $\ d$ denotes the 4 orbitals that we have. However, since we would like to calculate the Local Chern Marker (LCM), we would need to transform the above momentum space Hamiltonian into the real space. Using the fact that:

```math
c_{k, \alpha} = \frac{1}{N} \sum_{xy, \alpha} e^{-ik_x x} e^{-i k_y y}
```

where α denotes the different orbitals, and the relation that

```math
    \frac{1}{N} \sum_{k_x, k_y}e^{-ik_x(x-x')}e^{-ik_y(y-y')} = \delta_{xx'} \delta_{yy'}
```
we would then get the following $\ \textbf{real space}$ Hamiltonian

```math
    H_0 = \sum_{\mathbf{R}} [\sum_{s = a,b} (m+t_1) c^{\dagger}_{\mathbf{R}, s}c_{\mathbf{R}, s} + \sum_{p = c,d} (-m-t_1)c^{\dagger}_{\mathbf{R}, p}c_{\mathbf{R}, p}]\\
```

```math
 H_1 = \sum_{R} \sum_{\mu = \pm x, \pm y} \sum_{s = a,b; p = c,d} (\frac{-t_1}{4} c^{\dagger}_{\mathbf{R+l_\mu}, s} c_{\mathbf{R}, s}+\frac{t_1}{4}  c^{\dagger}_{\mathbf{R+l_\mu}, p} c_{\mathbf{R}, p} + e^{i\theta_\mu} \frac{t_2}{2} c^{\dagger}_{\mathbf{R+l_\mu}, s} c_{\mathbf{R}, p} + e^{-i \theta_\mu} \frac{t_2}{2}c^{\dagger}_{\mathbf{R+l_\mu}, (s;p)} c_{\mathbf{R}, (p;s)}) \\
```
```math
    H = H_0 + H_1
```
Some notations definitions:
* Here $\ s$ indicates the orbital $\ a \text{ and } b$ while $\ p$ indicates the orbital $\ c \text{ and } d$.
* $\ \mathbf{R}$ indicates the coordinate to one of the sites, in which $\ \mathbf{R} = (x,y)$.
* $\ \theta_\mu$ is the angle between the y-axis and the direction of the hopping $\ \mathbf{l_\mu}$ as also seen in Fig.(1).Since we only consider the nearest neighbour hopping then the value of $\ \theta_\mu$ would be $\ 0, \pm\pi/2, \pi$ only.
* $\ \frac{t_1}{4}$ indicates the hopping between the $\ s$-orbitals while $\ \frac{-t_1}{4}$ indicates the hopping between the $\ p$-orbitals.
* $\ \frac{t_2}{2}$ indicates the hopping between $\ s$ and $\ p$  orbitals and vice versa.
* $\ m$ indicates the "mass" of the orbital
* $\ H$ is the total real space Hamiltonian 

The above Hamiltonian produce a dispertion relation

```math
  E_k = \pm \sqrt{d_x^2+d_y^2+d_z^2}
```

With the following coefficient form

```math
    d_z^2=m^2 +t_1^2 \sin^4{(k_x/2)} + t_1^2 \sin^4{(k_y/2)}+2mt_1 \sin^2{(k_x/2)}+2mt_1 \sin^2{(k_y/2)}+t_1^2 \sin^2{(k_x/2)} \sin^2{(k_y/2)}
```
```math
  d_x^2 = t_2^2 \sin^2{k_x}
```
```math
  d_y^2 = t_2^2 \sin^2{k_y}
```

# LCM Calculation


# Localization Length Calculation
To actually check whether the phase transition of the LCM phase diagram indicates a topological phase transition, we attempt to relate it with the localization length or in other words how the disorder affects the quantum transport. So we start by a weird problem already: we have a 2D system, how do we calculate the localization since transfer matrices could only be calculated in 1D chain? This problem is answered in [Sarma's paper](https://iopscience.iop.org/article/10.1088/0022-3719/14/6/003) where we transform the problem into a quasi-1D problem which means that the system is approximately 1D but not actually 1D (weird right, but thats the creative part!). This happens for example if you see a cylinder with much much longer axial length compared to the radius, and if you see this very long radius from far away, what you see is just a 1D strip.

We define the quasi-1D localization length as $\ \rho_q1D$ while also defining the dimensionless quasi-1D localization length as $\ \Lambda = \rho_{q1D}/L$. By scaling analysis \textcolor{red}{(haven't understood what this is)} we then could get the localization of the original 2D shape system.
```math
  O_M = \prod_{n=1}^{M} T_n
```
According to Oseledec's theorem, the limit of
```math
  P = \lim_{M\to\infty} (O^{\dagger}_{M} O_{M})^{1/2M}
```
would have the eigenvalues of $\ \{e^{(\nu_1)}, e^{(-\nu_1)}, \ldots, e^{(\nu_s)}, e^{(-\nu_s)}\}$, in which $\ \nu_i \geq \nu_{i+1} \geq 0$ with $\ i=1,2, \ldots, s$
```math
    O_M = UR
```
```math
    \nu_i = \lim_{M\to\infty} \frac{\ln{(R)_{i,i}}}{M}
```




