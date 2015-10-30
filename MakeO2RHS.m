function C = MakeO2RHS(para, gas_o2, gas_no, u, v, r_coeff_o2, a, ind_r2,...
    ind_r3, ind_r4, nr_i_01, nr_i_12, nr_i_23)
  % Extra reaction terms for the RHS
  R_NO_23 = gas_no.R_max .* v(ind_r2 + 1 : ind_r3) ./...
      (v(ind_r2 + 1 : ind_r3) + gas_o2.Km_eNOS);
  app_Km_34 = para.Km .* (1 + u(ind_r3 + 1 : ind_r4) ./ gas_no.C_ref);
  R_O2_34 = gas_o2.Q_max_vw .* v(ind_r3 + 1 : ind_r4) ./...
      (v(ind_r3 + 1 : ind_r4) + app_Km_34);
  app_Km_45 = para.Km .* (1 + u(ind_r4 + 1 : end) ./ gas_no.C_ref);
  R_O2_45 = gas_o2.Q_max_t .* v(ind_r4 + 1 : end) ./...
      (v(ind_r4 + 1 : end) + app_Km_45);
  % RHS by layer
  C_01 = gas_o2.P * ones(nr_i_01, 1);
  C_12 = zeros(nr_i_12, 1);
  C_23 = r_coeff_o2 .* R_NO_23;
  C_34 = r_coeff_o2 .* R_O2_34;
  C_45 = r_coeff_o2 .* R_O2_45;
  % Overall RHS
  C = [gas_o2.P; C_12; C_23; C_34; C_45];
end
