clear;
clf;
u_ans = zeros(276, 5);
v_ans = zeros(276, 5);
tic();
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
  omega = 1.25; % factor for successive over relaxation method
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
  u_new = zeros(nr, 1);  % solution for NO
  v = [gas_o2.P * ones(nr_01, 1); zeros(nr - nr_01, 1)];  % solution for O2
  v_new = [gas_o2.P * ones(nr_01, 1); zeros(nr - nr_01, 1)];  % solution for O2
  r = linspace(0, para.R, nr);
  % Display interface positions
  fprintf('r0: %.1f\nr1: %.1f\nr2: %.1f\nr3: %.1f\nr4: %.1f\nr5: %.1f\n',...
      r(1), r(ind_r1), r(ind_r2), r(ind_r3), r(ind_r4), r(end));

  % Set up coefficient for R terms for clarity
  r_coeff_no = h * h / gas_no.d_coeff;
  r_coeff_o2 = h * h / gas_o2.d_coeff / para.alpha;

  %% LHS
  a = h / 2 ./ r (2 : end - 1);
  % Create LHS for NO
  diags_M_NO = [[1 - a'; 0; 0],...
                -2 / omega * ones(nr, 1)];
  diag_N_NO = [(2 - 2 / omega) * ones(nr, 1),...
               [0; -2; -1 - a']];
  M_NO = spdiags(diags_M_NO, [-1; 0], nr, nr);
  N_NO = spdiags(diag_N_NO, [0; 1], nr, nr);
  N_NO(end, end - 1) = -2;
  % Create LHS for O2
  diags_M_O2 = [[1 - a(ind_r1 : end)'; 0; 0],...
                [1; -2 / omega * ones(nr_i - nr_i_01 + 1, 1)]];
  diag_N_O2 = [[0; (2 - 2 / omega) * ones(nr - nr_i_01 - 1, 1)],...
               [0; 0; -1 - a(ind_r1 : end)']];
  M_O2 = spdiags(diags_M_O2, [-1; 0], nr - nr_i_01, nr - nr_i_01);
  N_O2 = spdiags(diag_N_O2, [0; 1],  nr - nr_i_01, nr - nr_i_01);
  N_O2(end, end - 1) = -2;
  % Solve for an initial O2 profile
  G = MakeO2RHS(para, gas_o2, gas_no, u, v, r_coeff_o2, a, ind_r2, ind_r3,...
      ind_r4, nr_i_01, nr_i_12, nr_i_23);
  v(ind_r1 : end) = M_O2 \ (N_O2 * v(ind_r1 : end) + G);

  tolerance = 1e-7;
  is_unsteady = true;
  % is_unsteady = false;
  while is_unsteady
    % Solving for NO
    F = MakeNORHS(para, gas_o2, gas_no, u, v, r_coeff_no, a, ind_r1, ind_r2,...
        ind_r3, ind_r4, nr_i_01, nr_i_12, nr_i_23, lambda_core);
    u_new = M_NO \ (N_NO * u + F);
    % Solving for O2
    G = MakeO2RHS(para, gas_o2, gas_no, u, v, r_coeff_o2, a, ind_r2, ind_r3,...
        ind_r4, nr_i_01, nr_i_12, nr_i_23);
    v_new(ind_r1 : end) = M_O2 \ (N_O2 * v(ind_r1 : end) + G);

    disp(max(abs(u_new - u)));
    % Checks if steady state is reached
    if max(abs(u_new - u)) < tolerance
      is_unsteady = false;
    end
    u = u_new;
    v = v_new;
  end
  u_ans(:, cfl_width) = u;
  v_ans(:, cfl_width) = v;
end
disp(toc());
% Write to data file for post processing
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
