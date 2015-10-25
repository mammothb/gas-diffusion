%===============================================================================
% Makes a n-by-n LHS matrix which solves for Dirichlet (left) BC and Neumann
% (right) BC. Can change between one-sided approximation or central
% approximation for the Neumann BC
% \param n size of one side of the matrix
% \param a the terms used to create the diagonals
% \return A the LHS matrix
%===============================================================================
function A = MakeDirichletNeumannMat(n, a)
  % one-sided
  diags = [[1 - a(2 : end)'; -4; 0], [-2 * ones(n - 1, 1); 3],...
      [0; 1 + a(1 : end)']];
  % imaginary node + correction
  % diags = [[1 - a(2 : end)'; 1; 0], [-2 * ones(n - 1, 1); -1],...
  %     [0; 1 + a(1 : end)']];
  A = spdiags(diags, [-1; 0; 1], n, n);
  A(end, end - 2) = 1;  % extra term for one-sided approximation of u'
end
