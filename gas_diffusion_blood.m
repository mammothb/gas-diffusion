clear;
clf;
u_ans = zeros(276, 5);
v_ans = zeros(276, 5);

for cfl_width = 1 : 5
  %% Model parameters
  para = Parameters();  % general parameters
  gas_no = Gas('NO', para);  % create object for NO gas
  gas_o2 = Gas('O2', para);  % create object for O2 gas
  cfl = cfl_width;  % [um], CFL width
  lambda_core = para.lambda_b / 2 * (1 + para.int_r * para.int_r /...
      (para.int_r - cfl) / (para.int_r - cfl));

  %% Simulation parameters
  h = 0.5;  % [um], space step
  if mod(para.R, h) > 1e-20
    error('Domain and space step incompatible');
  end
  nr = round(para.R / h) + 1;  % number of nodes in r direction
  nr_i = nr - 2;  % number of internal nodes
  ind_r1 = (para.int_r - cfl) / h + 1;  % index denoting end of RBC core
  ind_r2 = para.int_r / h + 1;  % index denoting end of vessel interior
  ind_r3 = ind_r2 + para.len_EC / h;  % index denoting end of EC layer
  ind_r4 = ind_r3 + para.len_VW / h;  % index denoting end of VW layer
  nr_01 = ind_r1;  % number of nodes for r0 < r < r1
  nr_12 = ind_r2 - ind_r1 + 1;  % number of nodes for r1 < r < r2 (include r = r1)
  nr_23 = ind_r3 - ind_r2 + 1;  % number of nodes for r2 < r < r3 (include r = r2)
  nr_34 = ind_r4 - ind_r3 + 1;  % number of nodes for r3 < r < r4 (include r = r3)
  nr_45 = nr - ind_r4 + 1;  % number of nodes for r4 < r < r5 (include r = r4)
  nr_i_01 = nr_01 - 1;  % size of matrix for r0 < r <= r1
  nr_i_12 = nr_12 - 1;  % size of matrix for r1 < r <= r2
  nr_i_23 = nr_23 - 1;  % size of matrix for r2 < r <= r3
  nr_i_34 = nr_34 - 1;  % size of matrix for r3 < r <= r4
  nr_i_45 = nr_45 - 1;  % size of matrix for r4 < r <= r5

  %% Initialization
  u = zeros(nr, 1);  % solution for NO
  u_new = zeros(nr, 1);  % auxillary variable, u_prime = u'
  v = [gas_o2.P * ones(nr_01, 1); zeros(nr - nr_01, 1)];  % solution for O2
  r = linspace(0, para.R, nr);
  fprintf('r0: %.1f\nr1: %.1f\nr2: %.1f\nr3: %.1f\nr4: %.1f\nr5: %.1f\n',...
      r(1), r(ind_r1), r(ind_r2), r(ind_r3), r(ind_r4), r(end));

  % Set up coefficient for R terms for clarity
  r_coeff_no = h * h / gas_no.d_coeff;
  r_coeff_o2 = h * h / gas_o2.d_coeff / para.alpha;

  %% LHS
  a = h / 2 ./ r (2 : end - 1);

  diags_O2 = [[zeros(nr_01 - 1, 1); 1 - a(ind_r1 : end)'; -4; 0],...
              [ones(nr_01, 1); -2 * ones(nr_i - nr_i_01, 1); 3],...
              [zeros(nr_01, 1); 0; 1 - a(ind_r1 : end)']];
  A_O2 = spdiags(diags_O2, [-1; 0; 1], nr, nr);
  A_O2(end, end - 2) = 1;

  C = MakeO2RHS(para, gas_o2, gas_no, u, v, r_coeff_o2, a, ind_r2, ind_r3,...
      ind_r4, nr_i_01, nr_i_12, nr_i_23);
  v = A_O2 \ C;
  tic();
  tol = 1e-6;
  flag = true;
  while flag
    for jj = 1 : nr
      rhs = GetNORHS(para, gas_no, gas_o2, ind_r1, ind_r2, ind_r3, ind_r4,...
          jj, u(jj), v(jj), lambda_core);
      rhs_o2 = GetO2RHS(para, gas_no, gas_o2, ind_r1, ind_r2, ind_r3, ind_r4,...
          jj, u(jj), v(jj), lambda_core);
      if jj == 1
        u_new(jj) = -0.5 * (-2 * u(jj + 1) + r_coeff_no * rhs);
      elseif jj == nr
        u_new(jj) = -0.5 * (-2 * u(jj - 1) + r_coeff_no * rhs);
      else
        u_new(jj) = -0.5 * (-(1 - a(jj - 1)) * u_new(jj - 1) -...
            (1 + a(jj - 1)) * u(jj + 1) + r_coeff_no * rhs);
      end
      if jj <= ind_r1
        v_new(jj) = gas_o2.P;
      elseif jj == nr
        v_new(jj) = -0.5 * (-2 * v(jj - 1) + r_coeff_o2 * rhs_o2);
      else
        v_new(jj) = -0.5 * (-(1 - a(jj - 1)) * v_new(jj - 1) -...
            (1 + a(jj - 1)) * v(jj + 1) + r_coeff_o2 * rhs_o2);
      end
    end  % jj
    disp(max(abs(u_new - u)));
    if max(abs(u_new - u)) < tol
      flag = false;
    end
    u = u_new;
    v = v_new;
    % C = MakeO2RHS(para, gas_o2, gas_no, u, v, r_coeff_o2, a, ind_r2, ind_r3,...
    %     ind_r4, nr_i_01, nr_i_12, nr_i_23);
    % v = A_O2 \ C;
  end
  disp(toc());
  u_ans(:, cfl_width) = u;
  v_ans(:, cfl_width) = v;
end
% subplot(2, 1, 1);
% plot(r(1 : ind_r1), u(1 : ind_r1),...
%     r(ind_r1 + 1 : ind_r2), u(ind_r1 + 1 : ind_r2),...
%     r(ind_r2 + 1 : ind_r3), u(ind_r2 + 1 : ind_r3),...
%     r(ind_r3 + 1 : ind_r4), u(ind_r3 + 1 : ind_r4),...
%     r(ind_r4 + 1 : end), u(ind_r4 + 1 : end));
% legend('RBC', 'CFL', 'EC', 'VW', 'T');
% subplot(2, 1, 2);
% plot(r(1 : ind_r1), v(1 : ind_r1),...
%     r(ind_r1 + 1 : ind_r2), v(ind_r1 + 1 : ind_r2),...
%     r(ind_r2 + 1 : ind_r3), v(ind_r2 + 1 : ind_r3),...
%     r(ind_r3 + 1 : ind_r4), v(ind_r3 + 1 : ind_r4),...
%     r(ind_r4 + 1 : end), v(ind_r4 + 1 : end));
dlmwrite('data.dat', [r', u_ans, v_ans], 'delimiter', ' ');
subplot(2, 1, 1);
plot(r, u_ans(:, 1),...
     r, u_ans(:, 2),...
     r, u_ans(:, 3),...
     r, u_ans(:, 4),...
     r, u_ans(:, 5));
legend('1', '2', '3', '4', '5');
subplot(2, 1, 2);
plot(r, v_ans(:, 1),...
     r, v_ans(:, 2),...
     r, v_ans(:, 3),...
     r, v_ans(:, 4),...
     r, v_ans(:, 5));
legend('1', '2', '3', '4', '5');
