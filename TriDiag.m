%===============================================================================
% Creates sparse tri-diagonal matrix of size n x n
%===============================================================================
function T = TriDiag(n, Tmain, Tsub, Tsup)
  T = spalloc(n, n, 3 * n);
  T_main = sparse(1 : n, 1 : n, Tmain * ones(1, n), n, n);
  T_sub = sparse(2 : n, 1 : n - 1, Tsub * ones(1, n - 1), n, n);
  T_sup = sparse(1 : n - 1, 2 : n, Tsup * ones(1, n - 1), n, n);
  T = T_main + T_sub + T_sup;
end
