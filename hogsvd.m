function [U, S, V, Tau, taumin, taumax, iso_classes] = hogsvd(A, m, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [U, S, V, Tau, taumin, taumax, m, iso_classes] = hogsvd(A, m, varargin)
%
% Computes the HO-GSVD of N=length(m) matrices with A=[A1; ... AN] and Ai
% of size m(i) x n. Returns U=[U1; ... UN], S=[S1; ... SN], and V s.t.
%       Ai=Ui*Si*V',
%       Ui=U(1+sum(m(1:i-1)):sum(m(1:i)),:),
%       Si=S(1+n*(i-1):n*i,:).
% If rank(A)<n, the matrix A is padded and the algorithm
% adds a matrix AN+1 and appends its dimension to m. Note that hogsvd(.)
% calls hocsd(.). See hocsd(.) for additional inputs and outputs.
%
% Optional keyword arguments (=default):
% RANK_TOL_A(=1e-14):       Rank tolerance for padding A
% ppi(=1):                  See hocsd(.)
% ZEROTOL(=1e-14):          See hocsd(.)
% EPS_REL_ISO(=1e-6):       See hocsd(.)
% DISABLE_WARNINGS(=false): Disable warnings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = size(A,2);
    N = length(m);
    if size(A,1)~=sum(m)
        error('size(A,1)=%d ~= sum(m)=%d.',size(A,1),sum(m));
    end
    args = parse_args(varargin{:});
    
    SA = svd(A);
    rank_A = sum(SA > args.RANK_TOL_A);
    rank_def_A = n - rank_A;
    if rank_def_A == 0
        Apad = A;
    else
        if ~args.DISABLE_WARNINGS
            warning('Provided rank-deficient A w. rank(A)=%d<n=%d. Padding A.', rank_A, n);
        end
        [~, ~, VA] = svd(A);
        Apad = [A; VA(:,end+1-rank_def_A:end)'];
        m = [m, rank_def_A]; % extend m
    end
    [Q, R] = qr(Apad, 0);
    [U, S, Z, Tau, taumin, taumax, iso_classes] = hocsd(Q, m, varargin{:});
    V = R'*Z;
    
    if rank_def_A > 0 % remove padding
        U(end+1-rank_def_A:end,:) = [];
        S(end-n+1:end,:) = [];
    end
end

function args = parse_args(varargin)
    check_POS = @(x) (x>0);
    check_EPS = @(x) ((x>0) && (x<1));
    p = inputParser;
    addParameter(p, 'RANK_TOL_A', 1e-14, check_EPS);
    addParameter(p, 'ppi', 1, check_POS);
    addParameter(p, 'ZEROTOL', 1e-14, check_EPS);
    addParameter(p, 'EPS_REL_ISO',	1e-6, check_EPS);
    addParameter(p, 'DISABLE_WARNINGS', false);
    parse(p, varargin{:});
    args = p.Results;
end