clear;
clf;
%% Parameters that do not change in the loop
% Model parameters
params = Parameters();  % general parameters
normal_cfl = (0.213 + 0.135) * params.int_r / 2.0;
rbc_core_radius = params.int_r - normal_cfl;
offset_radius = 1.0;
offset_angle = 1.0;
quarter_coordinates = GetQuarterCoordinates(offset_radius, offset_angle,...
    rbc_core_radius);
disp(quarter_coordinates);

%===============================================================================
% Simulation parameters/flags
%
% \param h [um], space step
% \param omega Factor for Successive over relaxation method
% \param tolerance Tolerance of relative error between iterations
% \param start_cfl Starting CFL width
% \param max_cfl [um], Maximum CFL width to test for
% \param which_scheme Determines which scheme to use for Neumann BCs
%        1: centered approx + correction term, first-order accurate
%        2: one-sided approx, second-order accurate
% \param cut_out_len Index deciding the amount of data to extract
% \param show_err Boolean for displaying relative error between iterations
% \param write_data Boolean to decide where results are written to data file
% \param write_small_data Boolean to decide whether we are going to cut out part
%        of the data for clearer analysis
% \param show_plot Boolean to decide if plot are to be drawn
%===============================================================================
h = 0.1;
omega = 1.9;
tolerance = 1e-6;
start_point = 1;
end_point = 4;
which_scheme = 2;
cut_out_len = 100.0;
show_err = false;
write_data = true;
write_small_data = true;
show_plot = true;

% Check if space step size is appropriate
if mod(params.R, h) > 1e-20, error('Domain and space step incompatible'); end

% Various compartments and domain space
nr = round(params.R / h) + 1;  % number of nodes in r direction
r = linspace(0, params.R, nr);

% Big matrix for all solutions
u_ans = zeros(nr, end_point - start_point + 1);
v_ans = zeros(nr, end_point - start_point + 1);

% Set up coefficient for R terms for clarity
r_coeff_no = h * h / params.no.d_coeff;
r_coeff_o2 = h * h / params.o2.d_coeff / params.alpha;

% a terms for LHS matrix
a = h / 2 ./ r(2 : end - 1);

% Create LHS for NO
[M_NO, N_NO] = MakeNOLHS(which_scheme, omega, a, nr);

tic();  % start stopwatch
for ii = start_point : end_point
  iteration = 0;
  cfl = params.int_r - quarter_coordinates(ii);  % [um], CFL width
  lambda_core = params.lambda_b / 2 * (1 + params.int_r * params.int_r /...
      (params.int_r - cfl) / (params.int_r - cfl));
  is_unsteady = true;  % while loop toggle

  %% Various compartments and domain space
  % Interface indexes
  ind_r1 = round((params.int_r - cfl) / h + 1);  % index for end of RBC core
  ind_r2 = round(params.int_r / h + 1);          % index for end of inner vessel
  ind_r3 = round(ind_r2 + params.len_EC / h);    % index for end of EC layer
  ind_r4 = round(ind_r3 + params.len_VW / h);    % index for end of VW layer
  % index denoting end of data to extract
  ind_extract = round(ind_r4 + cut_out_len / h);
  % Number of nodes in selected compartments
  nr_01 = round(ind_r1);                      % number of nodes for r0 < r <= r1
  nr_12 = round(ind_r2 - ind_r1);             % number of nodes for r1 < r <= r2
  nr_15 = round(nr - nr_01);                  % number of nodes for r1 < r <= r5
  % Index vectors for selected compartments for clarity
  r_01 = 1 : ind_r1;                      % r0 < r <= r1
  r_23 = ind_r2 + 1 : ind_r3;             % r2 < r <= r3
  r_34 = ind_r3 + 1 : ind_r4;             % r3 < r <= r4
  r_45 = ind_r4 + 1 : nr;                 % r4 < r <= r5

  %% Initialization
  u = zeros(nr, 1);                                         % solution for NO
  u_new = zeros(nr, 1);                                     % solution for NO
  v = [params.o2.P * ones(nr_01, 1); zeros(nr_15, 1)];      % solution for O2
  v_new = [params.o2.P * ones(nr_01, 1); zeros(nr_15, 1)];  % solution for O2

  % Display interface positions
  fprintf('CFL width: %d\n', cfl);
  fprintf('   r0    r1    r2    r3    r4    r5 \n');
  fprintf('%5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n',...
      r(1), r(ind_r1), r(ind_r2), r(ind_r3), r(ind_r4), r(end));

  % Create LHS for O2
  [M_O2, N_O2] = MakeO2LHS(which_scheme, omega, a, nr_15, ind_r1);
  % Solve for an initial O2 profile
  G = MakeO2RHS(params, u, v, r_coeff_o2, which_scheme, nr_12, r_23, r_34,...
      r_45);
  v(ind_r1 : end) = M_O2 \ (N_O2 * v(ind_r1 : end) + G);

  % Starts iterating to find the answer
  while is_unsteady
    iteration = iteration + 1;
    % Solving for NO
    F = MakeNORHS(params, u, v, r_coeff_no, which_scheme, r_01, nr_12, r_23,...
        r_34, r_45, lambda_core);
    u_new = M_NO \ (N_NO * u + F);
    % Solving for O2
    G = MakeO2RHS(params, u, v, r_coeff_o2, which_scheme, nr_12, r_23, r_34,...
        r_45);
    v_new(ind_r1 : end) = M_O2 \ (N_O2 * v(ind_r1 : end) + G);

    % Calculate relative errors
    u_relative_error = max(abs((u_new - u) ./ u_new));
    v_relative_error = max(abs((v_new - v) ./ v_new));
    % Output relative errors per iteration
    if show_err, fprintf('%.5d %.5d\n', u_relative_error, v_relative_error); end

    % Steady state check
    if u_relative_error < tolerance && v_relative_error < tolerance
      is_unsteady = false;
    end

    % Replace old solutions to prepare for next iteration
    u = u_new;
    v = v_new;
  end  % is_unsteady
  fprintf('Steady state reached in %d iterations\n', iteration);

  % Add solution to big solution matrix
  u_ans(:, ii) = u;
  v_ans(:, ii) = v;
end  % cfl_width
disp(toc());  % end stopwatch

% Write to data file for post processing
if write_data
  dlmwrite('data.dat', [r', u_ans, v_ans], 'delimiter', ' ');
  dlmwrite('data_small.dat', [r(1 : ind_extract)', u_ans(1 : ind_extract, :),...
      v_ans(1 : ind_extract, :)], 'delimiter', ' ');
end

% Plot all solutions
if show_plot
  subplot(2, 1, 1);
  plot(r, u_ans(:, 1),...
       r, u_ans(:, 2),...
       r, u_ans(:, 3),...
       r, u_ans(:, 4));
  legend('1', '2', '3', '4');
  subplot(2, 1, 2);
  plot(r, v_ans(:, 1),...
       r, v_ans(:, 2),...
       r, v_ans(:, 3),...
       r, v_ans(:, 4));
  legend('1', '2', '3', '4');
end
% fprintf('Peak: %d\nMean: %d\n', max(u_ans(:, 3)), mean(u_ans(:, 3)));
