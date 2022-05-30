function [U, S, Z, Tau, T, taumin, taumax] = hocsd(Q, N, m, n, ppi, ZEROTOL, EPS_REL_ISO, EPS_SVD_ISO)
% [U, S, Z, Tau, T, taumin, taumax] = hocsd(Q, N, m, n, ppi=1/N, ZEROTOL=1e-16, EPS_REL_ISO=1e-3, EPS_SVD_ISO=1e-3)
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
% dermined by selecting eigenvalues taui of T:=sum_i inv(Qi'*Qi+pi*I)/N
% that satisfy abs(taumax-taui)<=(taumax-taumin)*EPS_REL_ISO.
    enable_warnings = true;
    
    if ~exist('ppi','var') || isempty(ppi)
        ppi = 1/N;
    else
        if ppi<=0, error('Expected ppi>0, received ppi=%.4f', ppi), end
    end
    if (~exist('ZEROTOL','var') || isempty(ZEROTOL))
        ZEROTOL = 1e-14;
    end
    if ~exist('EPS_REL_ISO','var') || isempty(EPS_REL_ISO)
        EPS_REL_ISO = 1e-3;
    else
        if (EPS_REL_ISO<0 || EPS_REL_ISO>1), error('Expected EPS_REL_ISO=%.4f in [0,1]', EPS_REL_ISO), end
    end
    if (~exist('EPS_SVD_ISO','var')||isempty(EPS_SVD_ISO))
        EPS_SVD_ISO = 1e-3;
    else
        if (EPS_SVD_ISO<0 || EPS_SVD_ISO>1), error('Expected EPS_SVD_ISO=%.4f in [0,1]', EPS_SVD_ISO), end
    end
    if length(m)~=N, error('Expected length(m) == N'), end
    if size(Q, 1)~=sum(m), error('Expected size(Q, 1) == sum(m)'), end
    if size(Q, 2)~=n, error('Expected size(Q, 2) == n'), end
    if enable_warnings
        if norm(Q'*Q-eye(n))>=ZEROTOL, warning("Expected norm(Q'*Q-eye(n))/n=%e < ZEROTOL=%e",norm(Q'*Q-eye(n))/n,ZEROTOL), end
    end
    
    warn_eps_ortho_iso = 1e-6;
    warn_eps_eig_iso = 1e-6;
    warn_cond_Qi = 1e6;
    
    T = zeros(n, n);
    for i = 1 : N
        Qi = get_mat_from_stacked(Q, m, i);
        tmp = eye(n)/((Qi'*Qi) + eye(n,n)*ppi);
        T = T + tmp;
        if enable_warnings
            if cond(tmp)>=warn_cond_Qi, warning("For i=%d, cond(inv(Qi'Qi+piI))=%e\n",i,cond(tmp)), end
        end
    end

    T = (T+T') / 2 / N;
    [Z, Tau] = eig(T);
    [Tau, ind_Dt] = sort(diag(Tau), 'descend');
    Tau = diag(Tau);
    Z = Z(:, ind_Dt);

    taumin = 1/(1/N+ppi);
    taumax = (N-1)/N/ppi + 1/N/(1+ppi);
    ind_iso = abs(taumax*ones(n,1) - diag(Tau)) <= (taumax-taumin)*EPS_REL_ISO;
    if sum(ind_iso) > 0
        Z_iso = Z(:,ind_iso);
        n_iso = sum(ind_iso);
        Z_iso_rot = zeros(size(Z_iso));
        i_iso = 0;
        for j=1:N
            Qi = get_mat_from_stacked(Q, m, j);
            % TODO: Obtaining complex results (bad condition number)
            % [Xi,Si] = eig((Qi*Z_iso)'*(Qi*Z_iso));
            [~, Si, Xi] = svd(Qi*Z_iso); % using SVD instead
            ind_iso_i = diag(Si) >= 1-EPS_SVD_ISO;
            if sum(ind_iso_i)>0
                Z_iso_i = Z_iso*Xi(:,ind_iso_i);
                Z_iso_rot(:,1+i_iso:i_iso+sum(ind_iso_i)) = Z_iso_i;            
                i_iso = i_iso+sum(ind_iso_i);
                if (i_iso >= n_iso)
                    break;
                end
            end
            % TODO: Make sure no more isolated directions are in the
            % remainder of the loop.
        end
        Z(:,ind_iso) = Z_iso_rot;
        if (i_iso ~= n_iso), error('Found %d isolated direction through Qi, but expected %d from T', i_iso, n_iso), end % soft
        if enable_warnings
            if norm(eye(n)-Z'*Z)>warn_eps_ortho_iso, warning('Rotated Z is not orthogonal, norm(eye(n)-Z^T*Z)=%e.',norm(eye(n)-Z'*Z)), end
            if norm(T*Z-Z*Tau)>warn_eps_eig_iso, warning('Rotated is no EV matrix for T.'), end
        end
    end
    
    S = zeros(N*n, n);
    U = zeros(sum(m), n);
    for i = 1 : N
        Qi = get_mat_from_stacked(Q, m, i);
        Bi = Qi * Z;
        Si = diag(vecnorm(Bi, 2, 1));        
        for j=1:n
            if abs(Si(j,j)) > ZEROTOL
                U(1+sum(m(1:i-1)) : sum(m(1:i-1))+m(i), j) = Bi(:,j)/Si(j,j);
            else
                % TODO: Fill with orthogonal kernel basis
                if norm(Qi(:,j),2) > 0
                    U(1+sum(m(1:i-1)) : sum(m(1:i-1))+m(i), j) = ...
                        Qi(:,j)/norm(Qi(:,j),2);
                end
                Si(j,j) = 0;
            end
        end
        S(1+(i-1)*n : i*n, :) = Si;
    end
    
    % check
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
