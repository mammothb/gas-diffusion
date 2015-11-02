%===============================================================================
% This function creates the RHS vector required for solving for O2. The vector
% created is only from r1 <= r <= r5 since P_{O_2} is constant in the RBC
% compartment.
%
% \param para Struct containing various parameters for the model
% \param gas_o2 Struct containing various parameters related to O2 gas for the
%        model
% \param gas_no Struct containing various parameters related to NO gas for the
%        model
% \param u Vector containing values for C_{NO}
% \param v Vector containing values for P_{O_2}
% \param r_coeff_o2 Coefficient for the O2 reaction term to improve clarity
% \param which_scheme Scheme for handling Neumann BCs.
%        1: Centered approx + correction term, first-order accurate
%        2: One-sided approx, second-order accurate
% \param a Vector conatining values for 0.5 * h / r to improve clarity
% \param nr_12 Number of nodes for r1 < r <= r2
% \param r_23 Vector containing indexes for r2 < r <= r3
% \param r_34 Vector containing indexes for r3 < r <= r4
% \param r_45 Vector containing indexes for r4 < r <= r5
%===============================================================================
function C = MakeO2RHS(para, gas_o2, gas_no, u, v, r_coeff_o2, which_scheme,...
    nr_12, r_23, r_34, r_45)
  % Extra reaction terms for the RHS
  R_NO_23 = gas_no.R_max .* v(r_23) ./ (v(r_23) + gas_o2.Km_eNOS);
  app_Km_34 = para.Km .* (1 + u(r_34) ./ gas_no.C_ref);
  R_O2_34 = gas_o2.Q_max_vw .* v(r_34) ./ (v(r_34) + app_Km_34);
  app_Km_45 = para.Km .* (1 + u(r_45) ./ gas_no.C_ref);
  R_O2_45 = gas_o2.Q_max_t .* v(r_45) ./ (v(r_45) + app_Km_45);
  % RHS by layer
  C_12 = zeros(nr_12, 1);
  C_23 = r_coeff_o2 .* R_NO_23;
  C_34 = r_coeff_o2 .* R_O2_34;
  C_45 = r_coeff_o2 .* R_O2_45;
  % Overall RHS
  C = [gas_o2.P; C_12; C_23; C_34; C_45];
  % Modifies RHS based on selected scheme for Neumann BCs
  if which_scheme == 2, C(end) = 0; end
end
