function [U, S, Z, Tau, T, taumin, taumax] = hocsd(Q, N, m, n, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [U, S, Z, Tau, T, taumin, taumax] = hocsd(Q, N, m, n, varargin)
%
% Experimental version of the HO-CSD for rank-deficient matrices.
% Computes the HO-CSD of N matrices with Q = [Q1; Q2; ...; QN],
% where Qi is of size m(i) rows and n columns and satisfies
% Q'*Q=I.
%
% Returns U, S, Z such that Qi=Ui*Si*Z' with
% Ui=U(1+sum(m(1:i-1)):sum(m(1:i)),:) and Si=S(1+n*(i-1):n*i,:).
% 
% Tries to find an orthogonal basis for the isolated subspace, which is
% dermined by selecting eigenvalues taui of T:=sum_i inv(Qi'*Qi+ppi*I)/N
% that satisfy abs(taumax-taui)<=(taumax-taumin)*EPS_REL_ISO.
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
    if any(ind_iso)
        Z_iso = Z(:,ind_iso);
        n_iso = sum(ind_iso);
        Z_iso_new = zeros(size(Z_iso));
        i_iso = 0;
        for j=1:N
            Qi = get_mat_from_stacked(Q, m, j);
            [~, Si, Xi] = svd(Qi*Z_iso);
            % Handle different shapes of Si
            if (size(Si,2) == 1) || ((size(Si,1) == 1) && (size(Si,2) == 1))
                ind_iso_i = Si(1,1) >= 1-args.EPS_SVD_ISO;
            elseif size(Si,1) == 1
                ind_iso_i = Si(1,:) >= 1-args.EPS_SVD_ISO;
            else
                ind_iso_i = diag(Si) >= 1-args.EPS_SVD_ISO;
            end
            if any(ind_iso_i)
                Z_iso_new(:,1+i_iso:i_iso+sum(ind_iso_i)) = Z_iso*Xi(:,ind_iso_i);            
                i_iso = i_iso+sum(ind_iso_i);
                if (i_iso >= n_iso)
                    % TODO: Make sure no more isolated directions are in the
                    % remainder of the loop.
                    break;
                end
            end
        end
        Z(:,ind_iso) = Z_iso_new;
        if (i_iso ~= n_iso), error('Found %d isolated direction through Qi, but expected %d from T', i_iso, n_iso), end % soft
        if enable_warnings
            if norm(eye(n)-Z'*Z)>warn_eps_iso, warning('Rotated Z is not orthogonal, norm(eye(n)-Z^T*Z)=%e.',norm(eye(n)-Z'*Z)), end
            if norm(T*Z-Z*Tau)>warn_eps_iso, warning('Rotated is no EV matrix for T.'), end
        end
    end
    
    S = zeros(N*n, n);
    U = zeros(sum(m), n);
    for i = 1 : N
        Qi = get_mat_from_stacked(Q, m, i);
        Bi = Qi * Z;
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
    p = inputParser;
    addParameter(p, 'ppi', 1/N, check_TOL);
    addParameter(p, 'ZEROTOL', 1e-14, check_TOL);
    addParameter(p, 'EPS_REL_ISO', 1e-6, check_EPS);
    addParameter(p, 'EPS_SVD_ISO', 1e-6, check_EPS);
    parse(p, varargin{:});
    args = p.Results;
end