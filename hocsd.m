function [U, S, Z, Tau, T, taumin, taumax, iso_classes] = hocsd(Q, N, m, n, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [U, S, Z, Tau, T, taumin, taumax] = hocsd(Q, N, m, n, varargin)
%
% Experimental version of the HO-CSD for rank-deficient matrices.
% Computes the HO-CSD of N matrices with Q = [Q1; Q2; ...; QN],
% where Qi is of size m(i) rows and n columns and satisfies
% Q'*Q=I.
%
% Optional keyword arguments:
% ppi:          Scalar, default 1/N
% ZEROTOL:      Scalar, default 1e-14
% RANK_TOL_A:   Scalar, default 1e-14
% EPS_REL_ISO:  Scalar, default 1e-6
% EPS_SVD_ISO:  Scalar, default 1e-6
% GROUP_ISO:  	Boolean, true
% ACCELERATE:   Boolean, false
%
% Returns U, S, Z such that Qi=Ui*Si*Z' with
% Ui=U(1+sum(m(1:i-1)):sum(m(1:i)),:) and Si=S(1+n*(i-1):n*i,:).
% 
% Tries to find an orthogonal basis for the isolated subspace, which is
% dermined by selecting eigenvalues taui of T:=sum_i inv(Qi'*Qi+ppi*I)/N
% that satisfy abs(taumax-taui)<=(taumax-taumin)*EPS_REL_ISO. By setting
% the optional argument ACCELERATE to true, a simplified algorithm will
% determine the isolated directions using EPS_SVD_ISO, otherwise it will
% use a sequence of ever-thinner SVDs to determine the directions. The
% classes and indices corresponding to the isolated subspace are stored in
% the array iso_classes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    enable_warnings = true;
    if length(m)~=N, error('length(m)=%d ~= N=%d',length(m),N), end
    if sum(m)<n, error('sum(m)=%d < n=%d',sum(m),n), end
    if size(Q,1)~=sum(m), error('size(Q,1)=%d ~= sum(m)=%d',size(Q,1),sum(m)), end
    if size(Q,2)~=n, error('size(Q,2)=%d ~= n=%d',size(Q,2),n), end
    args = parse_args(N, varargin{:});
    if enable_warnings
        if norm(Q'*Q-eye(n))>=args.ZEROTOL, warning("Expected norm(Q'*Q-eye(n))/n=%e < ZEROTOL=%e",norm(Q'*Q-eye(n))/n,args.ZEROTOL), end
    end
    warn_eps_iso = 1e-6;
    warn_cond = 1e6;
    
    T = zeros(n, n);
    for i = 1 : N
        Qi = get_mat_from_stacked(Q, m, i);
        tmp = eye(n)/((Qi'*Qi) + eye(n,n)*args.ppi);
        T = T + tmp;
        if enable_warnings
            if cond(tmp)>=warn_cond, warning("For i=%d, cond(inv(Qi'Qi+piI))=%e\n",i,cond(tmp)), end
        end
    end

    T = (T+T') * (0.5/N);
    [Z, Tau] = eig(T);
    [Tau, ind_Dt] = sort(diag(Tau), 'descend');
    Tau = diag(Tau);
    Z = Z(:, ind_Dt);

    taumin = 1/(1/N+args.ppi);
    taumax = (N-1)/N/args.ppi + 1/N/(1+args.ppi);
    ind_iso = abs(taumax*ones(n,1) - diag(Tau)) <= (taumax-taumin)*args.EPS_REL_ISO;
    iso_classes = [];
    if any(ind_iso) % align Z_iso to standard RSVs of Qi
        Z_iso = Z(:,ind_iso);
        Z_iso_new = zeros(size(Z_iso));
        n_iso = sum(ind_iso);
        iso_classes = zeros(n_iso, 1);
        
        if args.ACCELERATE
            i_iso = 0;
            for j=1:N
                Qj = get_mat_from_stacked(Q, m, j);
                [~, Sj, Xj] = svd(Qj*Z_iso);
                if isvector(Sj)
                    Sj = Sj(1,1);
                elseif ~isscalar(Sj)
                    Sj = diag(Sj);
                end
                n_iso_j = sum(Sj > 1-args.EPS_SVD_ISO);
                if n_iso_j > 0
                    Z_iso_new(:,1+i_iso:i_iso+n_iso_j) = Z_iso*Xj(:,1:n_iso_j);
                    iso_classes(1+i_iso:i_iso+n_iso_j) = j;
                    i_iso = i_iso + n_iso_j;
                end
            end
            if i_iso ~= n_iso, error('Expected %d isolated directions, found %d. Increase EPS_SVD_ISO or disable ACCELERATE.', n_iso, i_iso); end
        else
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
        end
        
        Z(:,ind_iso) = Z_iso_new;
        if enable_warnings
            if norm(eye(n)-Z'*Z)>warn_eps_iso, warning('Rotated Z is not orthogonal, norm(eye(n)-Z^T*Z)=%e.',norm(eye(n)-Z'*Z,2)), end
            if norm(T*Z-Z*Tau)>warn_eps_iso, warning('Rotated Z is no EV matrix for T, norm(T*Z-Z*Tau)=%e. Decrease EPS_REL_ISO.',norm(T*Z-Z*Tau,2)), end
        end
    end
    if norm(eye(n)-Z'*Z) > warn_eps_iso; not_ortho=true; else; not_ortho=false; end
    
    S = zeros(N*n, n);
    U = zeros(sum(m), n);
    for i = 1 : N
        Qi = get_mat_from_stacked(Q, m, i);
        if not_ortho
            Bi = Qi / (Z');
        else
            Bi = Qi * Z;
        end
        Si = diag(vecnorm(Bi, 2, 1));        
        for j=1:n
            if Si(j,j) > args.ZEROTOL
                U(1+sum(m(1:i-1)):sum(m(1:i)), j) = Bi(:,j)*(1/Si(j,j));
            else
                % TODO: Fill with orthogonal kernel basis
                norm_Qi = norm(Qi(:,j),2);
                if norm_Qi > 0
                    U(1+sum(m(1:i-1)):sum(m(1:i)), j) = Qi(:,j)*(1/norm_Qi);
                end
            end
        end
        S(1+(i-1)*n:i*n, :) = Si;
    end
    
    if enable_warnings
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

function args = parse_args(N, varargin)
    check_TOL = @(x) (x>0);
    check_EPS = @(x) ((x>0) && (x<1));
    check = @(x) islogical(x);
    p = inputParser;
    addParameter(p, 'ppi',          1/N,	check_TOL);
    addParameter(p, 'ZEROTOL',      1e-14,	check_TOL);
    addParameter(p, 'EPS_REL_ISO',	1e-6,	check_EPS);
    addParameter(p, 'EPS_SVD_ISO',  1e-6,	check_EPS);
    addParameter(p, 'GROUP_ISO',    true,	check);
    addParameter(p, 'ACCELERATE',   false,	check);
    parse(p, varargin{:});
    args = p.Results;
end