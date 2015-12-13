%===============================================================================
% This function creates the LHS matrices required for solving O2
%
% \param which_scheme Determines which scheme to use for Neumann BCs
%        1: centered approx + correction term
%        2: one-sided approx
% \param omega Factor for SOR method
% \param a A term used in LHS matrices to improve clarity
% \param nr Number of nodes in r direction
% \param nr_15 Number of nodes between r1 and r5
% \param ind_r1 Index denoting end of RBC core
%===============================================================================
function [M, N] = MakeO2LHS(which_scheme, omega, a, nr, nr_15, ind_r1)
  switch which_scheme
  case 1
    D = spdiags([1; -2 * ones(nr_15, 1)], 0, nr_15 + 1, nr_15 + 1);
    L = spdiags([1 - a(ind_r1 : end)'; 2; 0], -1, nr_15 + 1, nr_15 + 1);
    U = spdiags([0; 0; 1 + a(ind_r1 : end)'], 1, nr_15 + 1, nr_15 + 1);
  case 2
    D = spdiags([1; -2 * ones(nr_15 - 1, 1); 3], 0, nr_15 + 1, nr_15 + 1);
    L = spdiags([1 - a(ind_r1 : end)'; -4; 0], -1, nr_15 + 1, nr_15 + 1);
    U = spdiags([0; 0; 1 + a(ind_r1 : end)'], 1, nr_15 + 1, nr_15 + 1);
    L(end, end - 2) = 1;
  otherwise
    error('Invalid Neumann BC scheme');
  end
  M = L + D ./ omega;
  N = (1 / omega - 1) .* D - U;
end
