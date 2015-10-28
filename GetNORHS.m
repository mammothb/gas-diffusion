function rhs = GetNORHS(para, gas_no, gas_o2, ind_r1, ind_r2, ind_r3, ind_r4,...
    index, u, v, lambda_core)
  if index <= ind_r1
    rhs = lambda_core * u;
  elseif index <= ind_r2
    rhs = 0;
  elseif index <= ind_r3
    rhs = -gas_no.R_max * v / (v + gas_o2.Km_eNOS);
  elseif index <= ind_r4
    rhs = para.lambda_vw * u;
  else
    rhs = para.lambda_t * u;
  end
end
