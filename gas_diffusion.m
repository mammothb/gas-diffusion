clear;
clf;
%% Parameters that do not change in the loop
% Model parameters
para = Parameters();  % general parameters
gas_no = Gas('NO', para);  % create object for NO gas
gas_o2 = Gas('O2', para);  % create object for O2 gas

% Simulation parameters
h = 0.5;  % [um], space step
omega = 1.75;  % factor for successive over relaxation method
tolerance = 1e-6;  % Tolerance for relative error for Gauss-Seidel
if mod(para.R, h) > 1e-20
  error('Domain and space step incompatible');
end

% Various compartments and domain space
nr = round(para.R / h) + 1;  % number of nodes in r direction
nr_i = nr - 2;  % number of internal nodes
r = linspace(0, para.R, nr);

% Big solution matrix
u_ans = zeros(nr, 5);
v_ans = zeros(nr, 5);

% Set up coefficient for R terms for clarity
r_coeff_no = h * h / gas_no.d_coeff;
r_coeff_o2 = h * h / gas_o2.d_coeff / para.alpha;

%% LHS
a = h / 2 ./ r(2 : end - 1);
% Create LHS for NO
D_NO = spdiags(-2 * ones(nr, 1), 0, nr, nr);
L_NO = spdiags([1 - a'; 2; 0], -1, nr, nr);
U_NO = spdiags([0; 2; 1 + a'], 1, nr, nr);
M_NO = L_NO + D_NO ./ omega;
N_NO = D_NO ./ omega - D_NO - U_NO;

tic();  % start stopwatch
for cfl_width = 1 : 5
  cfl = cfl_width;  % [um], CFL width
  lambda_core = para.lambda_b / 2 * (1 + para.int_r * para.int_r /...
      (para.int_r - cfl) / (para.int_r - cfl));

  is_unsteady = true;  % while loop toggle
  % is_unsteady = false;

  %% Various compartments and domain space
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
  % Display interface positions
  fprintf('CFL width: %d\n', cfl_width);
  fprintf('   r0    r1    r2    r3    r4    r5 \n');
  fprintf('%5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n',...
      r(1), r(ind_r1), r(ind_r2), r(ind_r3), r(ind_r4), r(end));

  % Create LHS for O2
  D_O2 = spdiags([1; -2 * ones(nr_i - nr_i_01 + 1, 1)], 0, nr - nr_i_01,...
      nr - nr_i_01);
  L_O2 = spdiags([1 - a(ind_r1 : end)'; 2; 0], -1, nr - nr_i_01, nr - nr_i_01);
  U_O2 = spdiags([0; 0; 1 + a(ind_r1 : end)'], 1, nr - nr_i_01, nr - nr_i_01);
  M_O2 = L_O2 + D_O2 ./ omega;
  N_O2 = D_O2 ./ omega - D_O2 - U_O2;
  % Solve for an initial O2 profile
  G = MakeO2RHS(para, gas_o2, gas_no, u, v, r_coeff_o2, a, ind_r2, ind_r3,...
      ind_r4, nr_i_01, nr_i_12, nr_i_23);
  v(ind_r1 : end) = M_O2 \ (N_O2 * v(ind_r1 : end) + G);
  % Starts iterating to find the answer
  while is_unsteady
    % Solving for NO
    F = MakeNORHS(para, gas_o2, gas_no, u, v, r_coeff_no, a, ind_r1, ind_r2,...
        ind_r3, ind_r4, nr_i_01, nr_i_12, nr_i_23, lambda_core);
    u_new = M_NO \ (N_NO * u + F);
    % Solving for O2
    G = MakeO2RHS(para, gas_o2, gas_no, u, v, r_coeff_o2, a, ind_r2, ind_r3,...
        ind_r4, nr_i_01, nr_i_12, nr_i_23);
    v_new(ind_r1 : end) = M_O2 \ (N_O2 * v(ind_r1 : end) + G);
    % Checks if steady state is reached
    u_relative_error = max(abs((u_new - u) ./ u_new));
    v_relative_error = max(abs((v_new - v) ./ v_new));
    % fprintf('%.5d %.5d\n', u_relative_error, v_relative_error);
    if u_relative_error < tolerance && v_relative_error < tolerance
      is_unsteady = false;
    end
    u = u_new;
    v = v_new;
  end  % is_unsteady
  % Add solution to big solution matrix
  u_ans(:, cfl_width) = u;
  v_ans(:, cfl_width) = v;
end  % cfl_width
disp(toc());  % end stopwatch

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
