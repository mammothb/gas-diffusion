%===============================================================================
% Creates sparse tri-diagonal matrix of the form
% | a b 0 ... 0 |
% | c a b ... 0 |
% | 0 c a ... 0 |
% | 0 ... c a b |
% | 0 ... 0 c a |
%===============================================================================
function T = TriDiag(n, a, b, c)
  T = spalloc(n, n, 3 * n);
  a_diag = sparse(1 : n, 1 : n, a * ones(1, n), n, n);
  b_diag = sparse(1 : n - 1, 2 : n, b * ones(1, n - 1), n, n);
  c_diag = sparse(2 : n, 1 : n - 1, c * ones(1, n - 1), n, n);
  T = a_diag + b_diag + c_diag;
end
