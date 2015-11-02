%===============================================================================
% This function creates the RHS vector required for solving for NO.
%
% \param para struct containing various parameters for the model
% \param gas_o2 struct containing various parameters related to O2 gas for the
%        model
% \param gas_no struct containing various parameters related to NO gas for the
%        model
% \param u vector containing values for C_{NO}
% \param v vector containing values for P_{O_2}
% \param r_coeff_no coefficient for the NO reaction term to improve clarity
% \param which_scheme Scheme for handling Neumann BCs.
%        1: Centered approx + correction term, first-order accurate
%        2: One-sided approx, second-order accurate
% \param a vector conatining values for 0.5 * h / r to improve clarity
% \param r_01 vector containing indexes for r0 < r <= r1
% \param nr_12 number of nodes for r1 < r <= r2
% \param r_23 vector containing indexes for r2 < r <= r3
% \param r_34 vector containing indexes for r3 < r <= r4
% \param r_45 vector containing indexes for r4 < r <= r5
% \param lambda_core lambda_core value needed to calculate NO reaction term in
%        r0 < r <=r1
%===============================================================================
function B = MakeNORHS(para, gas_o2, gas_no, u, v, r_coeff_no, which_scheme,...
    r_01, nr_12, r_23, r_34, r_45, lambda_core)
  % Extra reaction terms for the RHS
  R_NO_01 = lambda_core .* u(r_01);
  R_NO_23 = -gas_no.R_max .* v(r_23) ./ (v(r_23) + gas_o2.Km_eNOS);
  R_NO_34 = para.lambda_vw .* u(r_34);
  R_NO_45 = para.lambda_t .* u(r_45);
  % RHS by layer
  B_01 = r_coeff_no .* R_NO_01;
  B_12 = zeros(nr_12, 1);
  B_23 = r_coeff_no .* R_NO_23;
  B_34 = r_coeff_no .* R_NO_34;
  B_45 = r_coeff_no .* R_NO_45;
  % Overall RHS
  B = [B_01; B_12; B_23; B_34; B_45];
  switch which_scheme
  case 1
    B(1) = B(1) / 2;
  case 2
    B(1) = 0;
    B(end) = 0;
  otherwise
    error('Invalid Neumann BC scheme');
  end
end
