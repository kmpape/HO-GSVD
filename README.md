# HO-GSVD
Experimental Matlab implementation of the Higher-Order Generalized Singular Value (HO-GSVD) Decomposition. The implementation is based on the paper
"A Higher-Order Generalized Singular Value Decomposition for Rank Deficient Matrices" by Idris Kempf, Paul J. Goulart and Stephen R. Duncan. A pre-print
of the paper can be found at https://arxiv.org/abs/2102.09822.

The HO-GSVD is an extension of the standard SVD to $N\geq 2$ matrices. Given $N$ matrices $A_1,\dots,A_N$, the HO-GSVD decomposes each $A_i$ as $A_i=U_i\Sigma_i V^\text{T}$, $i = 1,\dots,N$, where $U_i\in\mathbb{R}^{m_i \times n}$ are the matrices of right basis vectors, 
the diagonal matrix $\Sigma_i\succeq 0$ contains the generalized singular values $\sigma_{i,k}$, 
and $V\in\mathbb{R}^{n\times n}$ with 
$\text{det}(V)\neq 0$ is the matrix of right basis vectors being shared among all factorizations. 

The HO-GSVD finds right-basis vectors $v_k$ that are exclusively used by one matrix $A_i$ 
($\sigma_{i,k}=1$ and $\sigma_{j,k}=0$ for $j\neq i$). 
If such a direction exist, the corresponding column of $U_i$, $i=1,\dots,N$, is orthogonal to all other columns. The HO-GSVD also finds right basis vectors $v_k$ that havel equal weight in each factorization, in the sense that 
$A_i^\text{T} A_i \tilde{v}_k = A_j^\text{T} A_j \tilde{v}_k$,
where $\tilde{v}_k = (A_i^\text{T} A_i + \pi I)^{-1} v_k = (A_j^\text{T} A_j + \pi I)^{-1} v_k$ 
and $\pi > 0$.
For more details, please see the corresponding paper.

If you use the HO-GSVD for published work or other projects, we encourage you to cite the accompanying paper.
```
@misc{hogsvd,
  doi       = {10.48550/ARXIV.2102.09822},  
  url       = {https://arxiv.org/abs/2102.09822},  
  author    = {Kempf, Idris and Goulart, Paul J. and Duncan, Stephen R.},  
  title     = {A Higher-Order Generalized Singular Value Decomposition for Rank Deficient Matrices},  
  publisher = {arXiv},  
  year      = {2021},
}
```
