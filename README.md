# RSCG_cpp
This package can calculate the elements of the Green's function:

```math
G_ij(σk) = ([σj I - A]^-1)_{ij},
```

with the use of the reduced-shifted conjugate gradient method
(See, Y. Nagai, Y. Shinohara, Y. Futamura, and T. Sakurai,[arXiv:1607.03992v2 or DOI:10.7566/JPSJ.86.014708]).
One can obtain ``G_{ij}(\sigma_k)`` with different frequencies ``\sigma_k``, simultaneously.

The matrix should be symmetric or hermitian.

This is to understand the method. 
I do not guarantee result with the use of this code. 