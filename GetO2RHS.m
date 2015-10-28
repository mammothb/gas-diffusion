function rhs = GetO2RHS(para, gas_no, gas_o2, ind_r1, ind_r2, ind_r3, ind_r4,...
    index, u, v, lambda_core)
  if index <= ind_r2
    rhs = 0;
  elseif index <= ind_r3
    rhs = gas_no.R_max * v / (v + gas_o2.Km_eNOS);
  elseif index <= ind_r4
    app_Km = para.Km * (1 + u / gas_no.C_ref);
    rhs = gas_o2.Q_max_vw * v / (v + app_Km);
  else
    app_Km = para.Km * (1 + u / gas_no.C_ref);
    rhs = gas_o2.Q_max_t * v / (v + app_Km);
  end
end
