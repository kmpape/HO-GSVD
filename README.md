# Higher-Order Generalized Singular Value Decomposition (HO-GSVD) and Higher-Order Cosine-Sine Decomposition (HO-CSD)
This is an experimental Matlab implementation of the Higher-Order Generalized Singular Value Decomposition based on the paper *A Higher-Order Generalized Singular Value Decomposition for Rank Deficient Matrices* by Idris Kempf, Paul J. Goulart and Stephen R. Duncan. A pre-print of the paper can be found at https://arxiv.org/abs/2102.09822.

The HO-GSVD is an extension of the standard SVD to $N\geq 2$ matrices. Given $N$ matrices $A_1,\dots,A_N$, the HO-GSVD decomposes each $A_i$ as $A_i=U_i\Sigma_i V^\text{T}$, $i = 1,\dots,N$, where $U_i\in\mathbb{R}^{m_i \times n}$ are the matrices of right basis vectors, 
the diagonal matrix $\Sigma_i\succeq 0$ contains the generalized singular values 
$\sigma_{i,k}$, 
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

# hogsvd.m
The HO-GSVD decomposes $N$ matrices $A_i\in\mathbb{R}^{m_i\times n}$, $i=1,\dots,N$, as $$A_i =U_i\Sigma_i V^T,$$ where the columns of $U_i\in\mathbb{R}^{m_i\times n}$ are referred to as *left basis vectors*, the diagonal elements of $\Sigma_i=\text{diag}(\sigma_{i,1},\dots,\sigma_{i,n})\in\mathbb{R}^{n\times n}$ as *generalized singular values* and the columns of $V\in\mathbb{R}^{n\times n}$ with $\text{det}(V)\neq 0$ as *right basis vectors*. Set $A:=[A_1^T, \dots, A_N^T]^T$ and $m:=[m_1,\dots,m_N]$, then call the HO-GSVD function as
```
[U, S, V, Tau, taumin, taumax, iso_classes] = hogsvd(A, m, varargin);
```
The function `hogsvd(A,m)` returns the *concatenated* factor matrices. To access the factor matrices of some $A_i$, use
```
Ui=U(1+sum(m(1:i-1)):sum(m(1:i)),:);
Si=S(1+n*(i-1):n*i,:);
```
The function `hogsvd(A,m)` also returns the diagonal matrix `Tau`, which contains the eigenvalues $\tau_k$, $k=1,\dots,n$, of $T_\pi$, where $$T_\pi:=\frac{1}{N}\sum_{i=1}^N(Q_i^T Q_i+\pi I)^{-1},$$
where the $Q_i$ are obtained from partitioning the left factor matrix of the thin QR factorization of $A$, i.e. from $QR=:A$ and $[Q_1^T,\dots,Q_N^T]^T:=Q$. Note that the present HO-GSVD implementation is based on the HO-CSD, i.e. `hogsvd(A,m)` calls `hocsd(Q,m)`.

Indices for which $\tau_k=\tau_\text{max}$, where $\tau_\text{max}$ is returned in `taumax`, are associated with the *isolated HO-GSVD* subspace, and the corresponding $\sigma_{i,k}=1$ for some $A_i$ and $\sigma_{j,k}=0$ for all other $A_j$, $j\neq i$. The matrix associated with the non-zero $\sigma_{i,k}$ is listed in `iso_classes`. Indices for which $\tau_k=\tau_\text{min}$, where $\tau_\text{min}$ is returned in `taumin`, are associated with the *common HO-GSVD* subspace, and the corresponding $\sigma_{i,k}=1/\sqrt{N}$ for all $A_i$. See the accompanying paper *A Higher-Order Generalized Singular Value Decomposition for Rank Deficient Matrices* for further details.

See folder `examples` for example usages.

# hocsd.m
The HO-CSD decomposes $N$ matrices $Q_i\in\mathbb{R}^{m_i\times n}$, $i=1,\dots,N$, satisfying $Q_1^TQ_1+\dots+Q_N^TQ_N=I,$ as $$Q_i =U_i\Sigma_i Z^T,$$ where $U_i\in\mathbb{R}^{m_i\times n}$, $\Sigma_i=\text{diag}(\sigma_{i,1},\dots,\sigma_{i,n})\in\mathbb{R}^{n\times n}$ and $Z\in\mathbb{R}^{n\times n}$ with $Z^T Z=I$. Set $Q:=[Q_1^T, \dots, Q_N^T]^T$ and $m:=[m_1,\dots,m_N]$, then call the HO-CSD function as
```
[U, S, Z, Tau, taumin, taumax, iso_classes] = hocsd(Q, m, varargin);
```
The function `hocsd(Q,m)` returns the *concatenated* factor matrices. To access the factor matrices of some $A_i$, use
```
Ui=U(1+sum(m(1:i-1)):sum(m(1:i)),:);
Si=S(1+n*(i-1):n*i,:);
```

See folder `examples` for example usages.
