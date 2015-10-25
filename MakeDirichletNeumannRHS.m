%===============================================================================
% Manipulates the RHS vector to add the appropriate boundary conditions, adds
% dirichlet BC to the first element in the vector, and neumann BC to the last
% element in the vector. Can change between one-sided approximation and Central
% approximation
% \param dirichlet_bc value of the Dirichlet BC
% \param neumann_bc value of the Neumann BC calculated using either one-sided or
%        central approximation
% \param a the term to be used with the Dirichlet BC
% \param h space step
% \param rhs the already populted RHS vector
% \return B the RHS vector after adjusting for Dirichlet and Neumann BC
%===============================================================================
function B = MakeDirichletNeumannRHS(dirichlet_bc, neumann_bc, a, h, rhs)
  B = rhs;
  B(1) = B(1) - (1 - a) * dirichlet_bc;
  % B(end) = B(end) - h * neumann_bc;
  B(end) = 2 * h * neumann_bc;
end
