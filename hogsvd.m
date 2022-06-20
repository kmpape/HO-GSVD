function [U, S, V, Q, R, Z, Tau, T, taumin, taumax, m, iso_classes] = hogsvd(A, N, m, n, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [U, S, V, Q, R, Z, Tau, T, taumin, taumax, m, iso_classes] = hogsvd(A, N, m, n, varargin)
%
% Computes the HO-GSVD of N matrices with A=[A1; A2; ...; AN] and Ai of
% size m(i) x n. If rank(A)<n, the matrix A is padded and the algorithm
% adds a matrix AN+1 and appends its dimension to m. 
% Returns U, S, V such that Ai=Ui*Si*V' with
% Ui=U(1+sum(m(1:i-1)):sum(m(1:i)),:) and Si=S(1+n*(i-1):n*i,:).
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
% Also returns the HO-CSD Qi=Ui*Si*Z', T and the eigenvalues of T, Tau, and
% the theoretical boundaries on the eigenvalues of T, taumin and taumax.
% The indices corresponding to the isolated subspace can be found in
% iso_classes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(m)~=N, error('length(m)=%d ~= N=%d',length(m),N), end
    if sum(m)<n, error('sum(m)=%d < n=%d',sum(m),n), end
    if size(A,1)~=sum(m), error('size(A,1)=%d ~= sum(m)=%d',size(A,1),sum(m)), end
    if size(A,2)~=n, error('size(A,2)=%d ~= n=%d',size(A,2),n), end
    args = parse_args(N, varargin{:});
    
    SA = svd(A);
    rank_A = sum(SA > args.RANK_TOL_A);
    rank_def_A = n - rank_A;
    if rank_def_A == 0
        Apad = A;
    else
        warning('Provided rank-deficient A w. rank(A)=%d<n=%d. Padding A.', rank_A, n);
        [~, ~, VA] = svd(A);
        Apad = [A; VA(:,end+1-rank_def_A:end)'];
        N = N + 1;
        m = [m, rank_def_A];
        % SApad = svd(Apad);
        % assert(sum(SApad >= args.RANK_TOL_A)==n);
    end
    [Q, R] = qr(Apad, 0);
%     [U, S, Z, Tau, T, taumin, taumax, iso_classes] = hocsd(Q, N, m, n,...
%         'ppi', args.ppi, 'ZEROTOL', args.ZEROTOL,...
%         'EPS_REL_ISO', args.EPS_REL_ISO, 'EPS_SVD_ISO', args.EPS_SVD_ISO);
    [U, S, Z, Tau, T, taumin, taumax, iso_classes] = hocsd(Q, N, m, n, varargin{:});
    V = R'*Z;
    
    if false && (rank_def_A>0) % remove padding
        U(end+1-rank_def_A:end,:) = [];
        S(end-n+1:end,:) = [];
    end
end

function args = parse_args(N, varargin)
    check_TOL = @(x) (x>0);
    check_EPS = @(x) ((x>0) && (x<1));
    check = @(x) islogical(x);
    p = inputParser;
    addParameter(p, 'RANK_TOL_A', 1e-14, check_TOL);
    % parameters for hocsd()
    addParameter(p, 'ppi',          1/N,	check_TOL);
    addParameter(p, 'ZEROTOL',      1e-14,	check_TOL);
    addParameter(p, 'EPS_REL_ISO',	1e-6,	check_EPS);
    addParameter(p, 'EPS_SVD_ISO',  1e-6,	check_EPS);
    addParameter(p, 'GROUP_ISO',    true,	check);
    addParameter(p, 'ACCELERATE',   false,	check);
    parse(p, varargin{:});
    args = p.Results;
end