clear;
clf;
%% Parameters that do not change in the loop
% Model parameters
para = Parameters();  % general parameters
gas_no = Gas('NO', para);  % create object for NO gas
gas_o2 = Gas('O2', para);  % create object for O2 gas

%===============================================================================
% Simulation parameters/flags
%
% \param h [um], space step
% \param omega Factor for Successive over relaxation method
% \param tolerance Tolerance of relative error between iterations
% \param max_cfl [um], Maximum CFL width to test for
% \param which_scheme Determines which scheme to use for Neumann BCs
%        1: centered approx + correction term, first-order accurate
%        2: one-sided approx, second-order accurate
% \param show_err Boolean for displaying relative error between iterations
% \param write_data Boolean to decide where results are written to data file
%===============================================================================
h = 0.5;
omega = 1.85;
tolerance = 1e-6;
max_cfl = 5;
which_scheme = 2;
cut_out_len = 100.0;
show_err = false;
write_data = true;
write_small_data = true;

% Check if space step size is appropriate
if mod(para.R, h) > 1e-20, error('Domain and space step incompatible'); end

% Various compartments and domain space
nr = round(para.R / h) + 1;  % number of nodes in r direction
r = linspace(0, para.R, nr);

% Big solution matrix
u_ans = zeros(nr, max_cfl);
v_ans = zeros(nr, max_cfl);

% Set up coefficient for R terms for clarity
r_coeff_no = h * h / gas_no.d_coeff;
r_coeff_o2 = h * h / gas_o2.d_coeff / para.alpha;

%% LHS
% Create LHS for NO
switch which_scheme
case 1
  D_NO = spdiags(-2 * ones(nr, 1), 0, nr, nr);
  L_NO = spdiags([ones(nr - 2, 1); 2; 0], -1, nr, nr);
  U_NO = spdiags([0; 2; ones(nr - 2, 1)], 1, nr, nr);
case 2
  D_NO = spdiags([-3; -2 * ones(nr - 2, 1); 3], 0, nr, nr);
  L_NO = spdiags([ones(nr - 2, 1); -4; 0], -1, nr, nr);
  U_NO = spdiags([0; 4; ones(nr - 2, 1)], 1, nr, nr);
  L_NO(end, end - 2) = 1;
  U_NO(1, 3) = -1;
otherwise
  error('Invalid Neumann BC scheme');
end
M_NO = L_NO + D_NO ./ omega;
N_NO = D_NO ./ omega - D_NO - U_NO;

tic();  % start stopwatch
for cfl_width = 1 : max_cfl
  cfl = cfl_width;  % [um], CFL width
  lambda_core = para.lambda_b / 2 * (1 + para.int_r * para.int_r /...
      (para.int_r - cfl) / (para.int_r - cfl));
  % while loop toggle
  is_unsteady = true;

  %% Various compartments and domain space
  % Interface indexes
  ind_r1 = (para.int_r - cfl) / h + 1;  % index denoting end of RBC core
  ind_r2 = para.int_r / h + 1;          % index denoting end of vessel interior
  ind_r3 = ind_r2 + para.len_EC / h;    % index denoting end of EC layer
  ind_r4 = ind_r3 + para.len_VW / h;    % index denoting end of VW layer
  % index denoting end of data to extract
  ind_extract = ind_r4 + cut_out_len / h;
  % Number of nodes in selected compartments
  nr_01 = ind_r1;                       % number of nodes for r0 < r <= r1
  nr_12 = ind_r2 - ind_r1;              % number of nodes for r1 < r <= r2
  nr_15 = nr - nr_01;                   % number of nodes for r1 < r <= r5
  % Index vectors for selected compartments for clarity
  r_01 = 1 : ind_r1;                    % r0 < r <= r1
  r_23 = ind_r2 + 1 : ind_r3;           % r2 < r <= r3
  r_34 = ind_r3 + 1 : ind_r4;           % r3 < r <= r4
  r_45 = ind_r4 + 1 : nr;               % r4 < r <= r5

  %% Initialization
  u = zeros(nr, 1);                                      % solution for NO
  u_new = zeros(nr, 1);                                  % solution for NO
  v = [gas_o2.P * ones(nr_01, 1); zeros(nr_15, 1)];      % solution for O2
  v_new = [gas_o2.P * ones(nr_01, 1); zeros(nr_15, 1)];  % solution for O2
  % Display interface positions
  fprintf('CFL width: %d\n', cfl_width);
  fprintf('   r0    r1    r2    r3    r4    r5 \n');
  fprintf('%5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n',...
      r(1), r(ind_r1), r(ind_r2), r(ind_r3), r(ind_r4), r(end));
  % Create LHS for O2
  switch which_scheme
  case 1
    D_O2 = spdiags([1; -2 * ones(nr_15, 1)], 0, nr_15 + 1, nr_15 + 1);
    L_O2 = spdiags([ones(nr_15, 1); 2; 0], -1, nr_15 + 1, nr_15 + 1);
    U_O2 = spdiags([0; 0; ones(nr_15, 1)], 1, nr_15 + 1, nr_15 + 1);
  case 2
    D_O2 = spdiags([1; -2 * ones(nr_15 - 1, 1); 3], 0, nr_15 + 1, nr_15 + 1);
    L_O2 = spdiags([ones(nr_15 - 1, 1); -4; 0], -1, nr_15 + 1, nr_15 + 1);
    U_O2 = spdiags([0; 0; ones(nr_15, 1)], 1, nr_15 + 1, nr_15 + 1);
    L_O2(end, end - 2) = 1;
  otherwise
    error('Invalid Neumann BC scheme');
  end
  M_O2 = L_O2 + D_O2 ./ omega;
  N_O2 = D_O2 ./ omega - D_O2 - U_O2;
  % Solve for an initial O2 profile
  G = MakeO2RHS(para, gas_o2, gas_no, u, v, r_coeff_o2, which_scheme, nr_12,...
      r_23, r_34, r_45);
  v(ind_r1 : end) = M_O2 \ (N_O2 * v(ind_r1 : end) + G);
  % Starts iterating to find the answer
  while is_unsteady
    % Solving for NO
    F = MakeNORHS(para, gas_o2, gas_no, u, v, r_coeff_no, which_scheme, r_01,...
        nr_12, r_23, r_34, r_45, lambda_core);
    u_new = M_NO \ (N_NO * u + F);
    % Solving for O2
    G = MakeO2RHS(para, gas_o2, gas_no, u, v, r_coeff_o2, which_scheme,...
        nr_12, r_23, r_34, r_45);
    v_new(ind_r1 : end) = M_O2 \ (N_O2 * v(ind_r1 : end) + G);
    % Calculate relative errors
    u_relative_error = max(abs((u_new - u) ./ u_new));
    v_relative_error = max(abs((v_new - v) ./ v_new));
    % Display relative errors
    if show_err, fprintf('%.5d %.5d\n', u_relative_error, v_relative_error); end
    % Steady state check
    if u_relative_error < tolerance && v_relative_error < tolerance
      is_unsteady = false;
    end
    % Replace old solutions to prepare for next iteration
    u = u_new;
    v = v_new;
  end  % is_unsteady
  % Add solution to big solution matrix
  u_ans(:, cfl_width) = u;
  v_ans(:, cfl_width) = v;
end  % cfl_width
disp(toc());  % end stopwatch
% Write to data file for post processing
if write_data
  dlmwrite('data.dat', [r', u_ans, v_ans], 'delimiter', ' ');
  dlmwrite('data_small.dat', [r(1 : ind_extract)', u_ans(1 : ind_extract, :),...
      v_ans(1 : ind_extract, :)], 'delimiter', ' ');
end
% Plot all solutions
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
