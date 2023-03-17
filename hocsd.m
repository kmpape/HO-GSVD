function [U, S, Z, Tau, taumin, taumax, iso_classes] = hocsd(Q, m, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [U, S, Z, Tau, taumin, taumax] = hocsd(Q, m, varargin)
%
% Experimental version of the HO-CSD for rank-deficient matrices.
% Computes the HO-CSD of N matrices with Q = [Q1; ... QN],
% where Qi is of size m(i) rows and n columns and satisfies
% Q'*Q=I. Returns U=[U1; ... UN], S=[S1; ... SN], and Z w. Z'Z=I such that
%       Qi=Ui*Si*Z',
%       Ui=U(1+sum(m(1:i-1)):sum(m(1:i)),:),
%       Si=S(1+n*(i-1):n*i,:).
% Also returns the eigenvalues Tau of T:=sum_i inv(Qi'*Qi+ppi*I)/N=Z*Tau*Z'
% in descending order. The dimension of the isolated subspace is
% n_iso=length(iso_classes), where iso_class(i) contains the index of the
% matrix Qi that has a unit generalised singular value.
%
% Optional keyword arguments (=default value):
% ppi(=1e-3):               Controls the condition number of Qi'*Qi+ppi*I
% ZEROTOL(=1e-14):          Tolerance for the geneneralised singular values
% EPS_REL_ISO(=1e-6):       Tolerance for the isolated subspace
% DISABLE_WARNINGS(=false): Disable warnings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WARN_EPS_ISO = 1e-6; % hard-coded
    WARN_COND = 1e6; % hard-coded
    N = length(m);
    n = size(Q,2);
    if sum(m)<n
        error('sum(m)=%d < n=%d. Rank(Q)=%d required.',sum(m),n,n);
    end
    if size(Q,1) ~= sum(m)
        error('size(Q,1)=%d ~= sum(m)=%d.',size(Q,1),sum(m));
    end
    
    args = parse_args(varargin{:});
    if ~args.DISABLE_WARNINGS
        if norm(Q'*Q-eye(n)) >= args.ZEROTOL
            warning("Expected norm(Q'*Q-eye(n))/n=%e < ZEROTOL=%e.",...
                norm(Q'*Q-eye(n))/n,args.ZEROTOL);
        end
    end
    
    % Compute eigenvectors and eigenvalues of T
    Rhat = zeros(n, N*n);
    sqrt_ppi = sqrt(args.ppi);
    for i = 1 : N
        Qi = get_mat_from_stacked(Q, m, i);
        [~, Rhati] = qr([Qi; eye(n)*sqrt_ppi], 0);
        Rhat(:, 1+(i-1)*n:i*n) = inv(Rhati);
        if ~args.DISABLE_WARNINGS
            if cond(Rhati)>=WARN_COND, warning("For i=%d, cond(Rhati)=%e\n",i,cond(Rhati)), end
        end
    end
    [Z, sqrt_Tau, ~] = svd(Rhat, 0);
    Tau = diag(diag(sqrt_Tau).^2/N);
    taumin = 1/(1/N+args.ppi); % theoretical min. of Tau
    taumax = (N-1)/N/args.ppi + 1/N/(1+args.ppi); % theoretical max. of Tau

    % Indices corresponding to the isolated subspace (eq. (6.6))
    ind_iso = abs(taumax*ones(n,1) - diag(Tau)) <= (taumax-taumin)*args.EPS_REL_ISO;
    
    % Align eigenvectors associated with the isolated subspace (Alg. 6.3)
    iso_classes = [];
    if any(ind_iso) % align Z_iso to standard RSVs of Qi
        Z_iso = Z(:,ind_iso);
        Z_iso_new = zeros(size(Z_iso));
        n_iso = sum(ind_iso);
        iso_classes = zeros(n_iso, 1);
        
        Z_iter = Z_iso;
        for i = 1:n_iso-1 % sort classes according to largest gain in the ith subspace
            all_S = zeros(N, 1);
            for j=1:N
                Qj = get_mat_from_stacked(Q, m, j);
                all_S(j) = norm(Qj*Z_iter, 2);
            end
            [~, ind_sorted] = sort(all_S, 'desc');
            iso_classes(i) = ind_sorted(1);
            Qiso = get_mat_from_stacked(Q, m, iso_classes(i));
            [~, ~, Xiso] = svd(Qiso*Z_iter);
            Z_iso_new(:,i) = Z_iter*Xiso(:,1);
            Z_iter = Z_iter*Xiso(:,2:end);
        end
        % Z_iter is a vector:
        all_S = zeros(N, 1);
        for j=1:N
            Qj = get_mat_from_stacked(Q, m, j);
            all_S(j,1) = norm(Qj*Z_iter,2);
        end
        [~, ind_sorted] = sort(all_S, 'desc');
        iso_classes(end) = ind_sorted(1);
        Z_iso_new(:,end) = Z_iter;
        
        Z(:,ind_iso) = Z_iso_new;
        if ~args.DISABLE_WARNINGS
            if norm(eye(n)-Z'*Z)>WARN_EPS_ISO
                warning('Rotated Z is not orthogonal, norm(eye(n)-Z^T*Z)=%e.',norm(eye(n)-Z'*Z,2))
            end
        end
    end
    if norm(eye(n)-Z'*Z) > WARN_EPS_ISO
        not_ortho=true;
    else
        not_ortho=false; 
    end
    
    S = zeros(N*n, n);
    U = zeros(sum(m), n);
    for i = 1 : N
        Qi = get_mat_from_stacked(Q, m, i);
        if not_ortho
            Bi = Qi / (Z');
        else
            Bi = Qi * Z;
        end
        Si = vecnorm(Bi, 2, 1);
        ind_pos = Si > args.ZEROTOL;
        U(1+sum(m(1:i-1)):sum(m(1:i)), ind_pos) = ...
            Bi(:,ind_pos) * diag(1./Si(ind_pos));
        
        % Process columns with zero GSVs (Lemma 2.5)
        nzero = sum(~ind_pos);
        if nzero > 0
            [UQi,SQi,~] = svd(Qi);
            ind_zero_i = diag(SQi)<=args.ZEROTOL;
            ni2 = sum(ind_zero_i);
            if ni2 == 0
                % Substitute normalised Qi columns
                Qitmp = Qi(:,~ind_pos);
                Qitmp_norm = vecnorm(Qitmp, 2, 1);
                Qitmp_norm(Qitmp_norm<=args.ZEROTOL) = 1;
                U(1+sum(m(1:i-1)):sum(m(1:i)), ~ind_pos) = Qitmp * diag(1./Qitmp_norm);
            else
                % Subsitute right SVs associated with zero SVs
                Ui2 = UQi(:,ind_zero_i);
                if ni2 < nzero
                    Ui2 = repmat(Ui2,[1,ceil(nzero/ni2)]);
                end
                U(1+sum(m(1:i-1)):sum(m(1:i)), ~ind_pos) = Ui2(:,1:nzero);
            end
        end
        S(1+(i-1)*n:i*n, :) = diag(Si);
    end
    
    if ~args.DISABLE_WARNINGS
        for i = 1 : N
            Ui = get_mat_from_stacked(U, m, i);
            Si = get_mat_from_stacked(S, repmat(n,1,N), i);
            Qi = get_mat_from_stacked(Q, m, i);
            if ~(sum(sum(abs(Ui*Si*Z'-Qi)))<1e-12)
                warning('HOCSD: reconstruction error=%e for matrix %d',...
                    sum(sum(abs(Ui*Si*Z'-Qi))), i);
            end
        end
    end
end

function [Ai, row_inds] = get_mat_from_stacked(A, m, i)
    assert(length(m) >= i);
    row_inds = 1+sum(m(1:i-1)) : sum(m(1:i));
    if ~exist('A','var')||isempty(A)
        Ai=[];
    else
        Ai = A(row_inds, :);
    end
end

function args = parse_args(varargin)
    check_POS = @(x) (x>=0);
    check_EPS = @(x) ((x>0) && (x<1));
    p = inputParser;
    addParameter(p, 'ppi',          1e-3,	check_POS);
    addParameter(p, 'ZEROTOL',      1e-14,	check_EPS);
    addParameter(p, 'EPS_REL_ISO',	1e-6,	check_EPS);
    addParameter(p, 'DISABLE_WARNINGS', false);
    parse(p, varargin{:});
    args = p.Results;
end