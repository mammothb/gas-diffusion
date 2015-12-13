%===============================================================================
% This function creates the LHS matrices required for solving NO
%
% \param which_scheme Determines which scheme to use for Neumann BCs
%        1: centered approx + correction term
%        2: one-sided approx
% \param omega Factor for SOR method
% \param a A term used in LHS matrices to improve clarity
% \param nr Number of nodes in r direction
%===============================================================================
function [M, N] = MakeNOLHS(which_scheme, omega, a, nr)
  switch which_scheme
  case 1
    D = spdiags(-2 * ones(nr, 1), 0, nr, nr);
    L = spdiags([1 - a'; 2; 0], -1, nr, nr);
    U = spdiags([0; 2; 1 + a'], 1, nr, nr);
  case 2
    D = spdiags([-3; -2 * ones(nr - 2, 1); 3], 0, nr, nr);
    L = spdiags([1 - a'; -4; 0], -1, nr, nr);
    U = spdiags([0; 4; 1 + a'], 1, nr, nr);
    L(end, end - 2) = 1;
    U(1, 3) = -1;
  otherwise
    error('Invalid Neumann BC scheme');
  end
  M = L + D ./ omega;
  N = (1 / omega - 1) .* D - U;
end
