function [U, S, V, Tau, T, taumin, taumax] = hogsvd(A, N, m, n, ppi, RANK_TOL_A,...
                                    ZEROTOL, eps_rel_iso, eps_svd_iso)
% [U, S, V, Tau, T, taumin, taumax] = hogsvd(A, N, m, n, ppi, RANK_TOL_A, ZEROTOL, eps_rel_iso, eps_svd_iso)
%
% Computes the HO-GSVD of N matrices with A=[A1; A2; ...; AN] and Ai of
% size m(i) x n. If rank(A)<n, the matrix A is padded. The double RANK_TOL_A
% is used to determine the rank of A.
%
% This function is based on the HO-CSD implementation. For the optional
% arguments ZEROTOL, eps_rel_iso, eps_svd_iso see hocsd(...).
%
% Returns U, S, Vs such that Ai=Ui*Si*V' with
% Ui=U(1+sum(m(1:i-1)):sum(m(1:i)),:) and Si=S(1+n*(i-1):n*i,:).
%
% Also returns T from the HO-CSD, and the eigenvalues of T, Tau, and the
% theoretical boundaries on the eigenvalues of T, taumin and taumax.
    if ~exist('RANK_TOL_A','var')
        RANK_TOL_A = 1e-14;
    else
        if RANK_TOL_A<0, error('Expected RANK_TOL_A>=0, received RANK_TOL_A=%.4f', RANK_TOL_A), end
    end
    if ~exist('ppi','var')
        ppi = [];
    else
        if ppi<=0, error('Expected ppi>0, received ppi=%.4f', ppi), end
    end
    if (~exist('ZEROTOL','var') || isempty(ZEROTOL))
        ZEROTOL = [];
    end
    if ~exist('eps_rel_iso','var')
        eps_rel_iso = [];
    else
        if (eps_rel_iso<0 || eps_rel_iso>1), error('Expected eps_rel_iso=%.4f in [0,1]', eps_rel_iso), end
    end
    if ~exist('eps_svd_iso','var')
        eps_svd_iso = [];
    else
        if (eps_svd_iso<0 || eps_svd_iso>1), error('Expected eps_svd_iso=%.4f in [0,1]', eps_svd_iso), end
    end
    if length(m)~=N, error('Expected length(m) == N'), end
    if size(A, 1)~=sum(m), error('Expected size(Aall, 1) == sum(m)'), end
    if size(A, 2)~=n, error('Expected size(Aall, 2) == n'), end
    
    SA = svd(A);
    rank_A = sum(SA > RANK_TOL_A);
    rank_def_A = n - rank_A;
    if rank_def_A == 0
        Apad = A;
    else
        warning('Provided rank-deficient A w. rank(A)=%d<n=%d. Padding A.', rank_A, n);
        [~, ~, VA] = svd(A);
        Apad = [A; VA(:,end+1-rank_def_A:end)'];
        N = N + 1;
        m =[m, rank_def_A];
        SApad = svd(Apad);
        assert(sum(SApad >= RANK_TOL_A)==n);
    end
    [Q, R] = qr(Apad, 0);
    [U, S, Z, Tau, T, taumin, taumax] = hocsd(Q, N, m, n, ppi, ZEROTOL, eps_rel_iso, eps_svd_iso);
    V = R'*Z;
    
    if false && (rank_def_A>0)
        U(end+1-rank_def_A:end,:) = [];
        S(end-n+1:end,:) = [];
    end
end
