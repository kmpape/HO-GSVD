function [Ai, row_inds] = get_mat_from_stacked(A, m, i)
% [Ai, row_inds] = get_mat_from_stacked(A, m, i)
%
% Given A=[A1; A2; ...; AN] with Ak of size size m(i) x n, this function
% returns Ai with the corresponding row indices row_inds.
    assert(length(m) >= i);
    rows_before = sum(m(1:i-1));
    row_inds = rows_before + 1 : rows_before + m(i);
    if ~exist('A','var')||isempty(A)
        Ai=[];
    else
        Ai = A(row_inds, :);
    end
end