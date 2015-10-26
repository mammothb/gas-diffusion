function B = MakeNORHS(para, gas_o2, gas_no, u, v, r_coeff_no, a, ind_r1,...
    ind_r2, ind_r3, ind_r4, nr_i_01, nr_i_12, nr_i_23, lambda_core)
  % Extra reaction terms for the RHS
  R_NO_01 = lambda_core .* u(2 : ind_r1);
  R_NO_23 = -gas_no.R_max .* v(ind_r2 + 1 : ind_r3) ./...
      (v(ind_r2 + 1 : ind_r3) + gas_o2.Km_eNOS);
  R_NO_34 = para.lambda_vw .* u(ind_r3 + 1 : ind_r4);
  R_NO_45 = para.lambda_t .* u(ind_r4 + 1 : end - 1);
  % RHS by layer
  B_01 = r_coeff_no .* R_NO_01;
  B_12 = zeros(nr_i_12, 1);
  B_23 = r_coeff_no .* R_NO_23;
  B_34 = r_coeff_no .* R_NO_34;
  B_45 = r_coeff_no .* R_NO_45;
  % Overall RHS
  B = [0; B_01; B_12; B_23; B_34; B_45; 0];
  % B(1) = B(1) - (1 - a(1)) * u(1);
  % B(end) = B(end) - (1 + a(end)) * u(end);
end
