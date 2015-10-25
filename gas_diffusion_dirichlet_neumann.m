clear;
clf;
%% Model parameters
para = Parameters();  % general parameters
gas_no = Gas('NO', para);  % create object for NO gas
gas_o2 = Gas('O2', para);  % create object for O2 gas
cfl = 2;  % [um], CFL width
lambda_core = para.lambda_b / 2 * (1 + para.int_r * para.int_r /...
    (para.int_r - cfl) / (para.int_r - cfl));

%% Simulation parameters
h = 0.005;  % [um], space step
if mod(para.R, h) > 1e-20
  error('Domain and space step incompatible');
end
nr = round(para.R / h) + 1;  % number of nodes in r direction
ind_r1 = (para.int_r - cfl) / h + 1;  % index denoting end of RBC core
ind_r2 = para.int_r / h + 1;  % index denoting end of vessel interior
ind_r3 = ind_r2 + para.len_EC / h;  % index denoting end of EC layer
ind_r4 = ind_r3 + para.len_VW / h;  % index denoting end of VW layer
nr_01 = ind_r1;  % number of nodes for r0 < r < r1
nr_12 = ind_r2 - ind_r1 + 1;  % number of nodes for r1 < r < r2 (include r = r1)
nr_23 = ind_r3 - ind_r2 + 1;  % number of nodes for r2 < r < r3 (include r = r2)
nr_34 = ind_r4 - ind_r3 + 1;  % number of nodes for r3 < r < r4 (include r = r3)
nr_45 = nr - ind_r4 + 1;  % number of nodes for r4 < r < r5 (include r = r4)
nr_i_01 = nr_01 - 1;  % size of matrix for r0 < r < r1
nr_i_12 = nr_12 - 1;  % size of matrix for r1 < r < r2
nr_i_23 = nr_23 - 1;  % size of matrix for r2 < r < r3
nr_i_34 = nr_34 - 1;  % size of matrix for r3 < r < r4
nr_i_45 = nr_45 - 2;  % size of matrix for r3 < r < r4

%% Initialization
u = zeros(nr, 1);  % solution for NO
v = [gas_o2.P * ones(nr_01, 1); zeros(nr - nr_01, 1)];  % solution for O2
% v = gas_o2.P * ones(nr, 1);  % solution for O2
r = linspace(0, para.R, nr);
r_01 = linspace(r(1), r(ind_r1), nr_01);
r_12 = linspace(r(ind_r1), r(ind_r2), nr_12);
r_23 = linspace(r(ind_r2), r(ind_r3), nr_23);
r_34 = linspace(r(ind_r3), r(ind_r4), nr_34);
r_45 = linspace(r(ind_r4), r(end), nr_45);
fprintf('r0: %.1f\nr1: %.1f\nr2: %.1f\nr3: %.1f\nr4: %.1f\nr5: %.1f\n',...
    r(1), r(ind_r1), r(ind_r2), r(ind_r3), r(ind_r4), r(end));
%===============================================================================
% All R_species terms have been treated as if they are shifted to the RHS
%===============================================================================
% set up coefficient for R terms for clarity
r_coeff_no = h * h / gas_no.d_coeff;
r_coeff_o2 = h * h / gas_o2.d_coeff / para.alpha;

%% Setting up all LHS matrix since they remain the same during the iterations
% RBC r0 < r < r1
a_01 = h / 2 ./ r(2 : ind_r1 - 1);
A_01 = MakeDirichletNeumannMat(nr_i_01, a_01);
% CFL r1 < r < r2
a_12 = h / 2 ./ r(ind_r1 + 1 : ind_r2 - 1);
A_12 = MakeDirichletNeumannMat(nr_i_12, a_12);
% EC r2 < r < r3
a_23 = h / 2 ./ r(ind_r2 + 1 : ind_r3 - 1);
A_23 = MakeDirichletNeumannMat(nr_i_23, a_23);
% VW r3 < r < r4
a_34 = h / 2 ./ r(ind_r3 + 1 : ind_r4 - 1);
A_34 = MakeDirichletNeumannMat(nr_i_34, a_34);
% T r4 < r < r5
a_45 = h / 2 ./ r(ind_r4 + 1 : end - 1);
diags_45 = [[1 - a_45(2 : end)'; 0], [-2 * ones(nr_i_45, 1)], [0;...
    1 + a_45(1 : end - 1)']];
A_45 = spdiags(diags_45, [-1; 0; 1], nr_i_45, nr_i_45);

%% Solving O2 to get initial partial pressure profile
% Neumann BC for continous mass flux using one-sided approx
% phi_v_12 = ForwardOneSideApprox(h, v(ind_r2), v(ind_r2 + 1), v(ind_r2 + 2));
phi_v_12 = CentralApprox(h, v(ind_r2 + 1), v(ind_r2 - 1));
% RHS
rhs_O2_12 = zeros(nr_i_12, 1);
C_12 = MakeDirichletNeumannRHS(v(ind_r1), phi_v_12, a_12(1), h, rhs_O2_12);
v(ind_r1 + 1 : ind_r2) = A_12 \ C_12;  % Solving
% Neumann BC for continuous mass flux using one-sided approx
% gamma_v_23 = ForwardOneSideApprox(h, v(ind_r3), v(ind_r3 + 1), v(ind_r3 + 2));
gamma_v_23 = CentralApprox(h, v(ind_r3 + 1), v(ind_r3 - 1));
% RHS
R_NO_23 = -gas_no.R_max .* v(ind_r2 + 1 : ind_r3) ./...
    (v(ind_r2 + 1 : ind_r3) + gas_o2.Km_eNOS);
rhs_O2_23 = r_coeff_o2 .* -R_NO_23;
C_23 = MakeDirichletNeumannRHS(v(ind_r2), gamma_v_23, a_23(1), h, rhs_O2_23);
v(ind_r2 + 1 : ind_r3) = A_23 \ C_23;  % Solving
% Neumann BC for continous mass flux using one-sided approx
% delta_v_34 = ForwardOneSideApprox(h, v(ind_r4), v(ind_r4 + 1), v(ind_r4 + 2));
delta_v_34 = CentralApprox(h, v(ind_r4 + 1), v(ind_r4 - 1));
% RHS
app_Km_34 = para.Km * (1 + u(ind_r3 + 1 : ind_r4) ./ gas_no.C_ref);
R_O2_34 = gas_o2.Q_max_vw .* v(ind_r3 + 1 : ind_r4) ./...
    (v(ind_r3 + 1 : ind_r4) + app_Km_34);
rhs_O2_34 = r_coeff_o2 .* R_O2_34;
C_34 = MakeDirichletNeumannRHS(v(ind_r3), delta_v_34, a_34(1), h, rhs_O2_34);
v(ind_r3 + 1 : ind_r4) = A_34 \ C_34;  % Solving
v(end) = v(end - 1);  % Neumann BC for zero mass flux
% RHS
app_Km_45 = para.Km .* (1 + u(ind_r4 + 1 : end - 1) ./ gas_no.C_ref);
R_O2_45 = gas_o2.Q_max_t .* v(ind_r4 + 1 : end - 1) ./...
    (v(ind_r4 + 1 : end - 1) + app_Km_45);
C_45 = r_coeff_o2 .* R_O2_45;
C_45(1) = C_45(1) - (1 - a_45(1)) * v(ind_r4);
C_45(end) = C_45(end) - (1 + a_45(end)) * v(end);
v(ind_r4 + 1 : end - 1) = A_45 \ C_45;  % Solving

for ii = 1 : 2
  u_old = u;  % copy previous solution to check for steady state
  v_old = v;
  %=============================================================================
  % RBC r0 < r < r1
  %=============================================================================
  % NO
  % Neumann BC for zero mass flux
  u(1) = u(2);
  % Neumann BC for continuous mass flux using one-sided approx
  % sigma_u_01 = ForwardOneSideApprox(h, u(ind_r1), u(ind_r1 + 1), u(ind_r1 + 2));
  % imaginary node
  sigma_u_01 = CentralApprox(h, u(ind_r1 + 1), u(ind_r1 - 1));
  % RHS
  R_NO_01 = lambda_core .* u(2 : ind_r1);
  rhs_no_01 = r_coeff_no .* R_NO_01;
  B_01 = MakeDirichletNeumannRHS(u(1), sigma_u_01, a_01(1), h, rhs_no_01);
  u(2:ind_r1) = A_01 \ B_01;  % Solving
  % O2
  % Nothing to be done for O2 since P_o2 is constant in this region

  %=============================================================================
  % CFL r1 < r < r2
  %=============================================================================
  % NO
  % Neumann BC for continuous mass flux using one-sided approx
  % phi_u_12 = ForwardOneSideApprox(h, u(ind_r2), u(ind_r2 + 1), u(ind_r2 + 2));
  phi_u_12 = CentralApprox(h, u(ind_r2 + 1), u(ind_r2 - 1));
  % RHS
  rhs_no_12 = zeros(nr_i_12, 1);
  B_12 = MakeDirichletNeumannRHS(u(ind_r1), phi_u_12, a_12(1), h, rhs_no_12);
  u(ind_r1 + 1 : ind_r2) = A_12 \ B_12;  % Solving
  % O2
  % Neumann BC for continous mass flux using one-sided approx
  % phi_v_12 = ForwardOneSideApprox(h, v(ind_r2), v(ind_r2 + 1), v(ind_r2 + 2));
  phi_v_12 = CentralApprox(h, v(ind_r2 + 1), v(ind_r2 - 1));
  % RHS
  rhs_o2_12 = zeros(nr_i_12, 1);
  C_12 = MakeDirichletNeumannRHS(v(ind_r1), phi_v_12, a_12(1), h, rhs_o2_12);
  v(ind_r1 + 1 : ind_r2) = A_12 \ C_12;  % Solving

  %=============================================================================
  % EC r2 < r < r3
  %=============================================================================
  % NO
  % Neumann BC for continuous mass flux using one-sided approx
  % gamma_u_23 = ForwardOneSideApprox(h, u(ind_r3), u(ind_r3 + 1), u(ind_r3 + 2));
  gamma_u_23 = CentralApprox(h, u(ind_r3 + 1), u(ind_r3 - 1));
  % NO production term
  R_NO_23 = -gas_no.R_max .* v(ind_r2 + 1 : ind_r3) ./...
      (v(ind_r2 + 1 : ind_r3) + gas_o2.Km_eNOS);
  % RHS
  rhs_no_23 = r_coeff_no .* R_NO_23;
  B_23 = MakeDirichletNeumannRHS(u(ind_r2), gamma_u_23, a_23(1), h, rhs_no_23);
  u(ind_r2 + 1 : ind_r3) = A_23 \ B_23;  % Solving
  % O2
  % Neumann BC for continuous mass flux using one-sided approx
  % gamma_v_23 = ForwardOneSideApprox(h, v(ind_r3), v(ind_r3 + 1), v(ind_r3 + 2));
  gamma_v_23 = CentralApprox(h, v(ind_r3 + 1), v(ind_r3 - 1));
  % RHS
  rhs_o2_23 = r_coeff_o2 .* -R_NO_23;
  C_23 = MakeDirichletNeumannRHS(v(ind_r2), gamma_v_23, a_23(1), h, rhs_o2_23);
  % Solving
  v(ind_r2 + 1 : ind_r3) = A_23 \ C_23;

  %=============================================================================
  % VW r3 < r < r4
  %=============================================================================
  % NO
  % Neumann BC for continous mass flux using one-sided approx
  % delta_u_34 = ForwardOneSideApprox(h, u(ind_r4), u(ind_r4 + 1), u(ind_r4 + 2));
  delta_u_34 = CentralApprox(h, u(ind_r4 + 1), u(ind_r4 - 1));
  R_NO_34 = para.lambda_vw .* u(ind_r3 + 1 : ind_r4);
  % RHS
  rhs_no_34 = r_coeff_no .* R_NO_34;
  B_34 = MakeDirichletNeumannRHS(u(ind_r3), delta_u_34, a_34(1), h, rhs_no_34);
  u(ind_r3 + 1 : ind_r4) = A_34 \ B_34;  % Solving
  % O2
  % Neumann BC for continous mass flux using one-sided approx
  % delta_v_34 = ForwardOneSideApprox(h, v(ind_r4), v(ind_r4 + 1), v(ind_r4 + 2));
  delta_v_34 = CentralApprox(h, v(ind_r4 + 1), v(ind_r4 - 1));
  app_Km_34 = para.Km * (1 + u(ind_r3 + 1 : ind_r4) ./ gas_no.C_ref);
  R_O2_34 = gas_o2.Q_max_vw .* v(ind_r3 + 1 : ind_r4) ./...
      (v(ind_r3 + 1 : ind_r4) + app_Km_34);
  % RHS
  rhs_o2_34 = r_coeff_o2 .* R_O2_34;
  C_34 = MakeDirichletNeumannRHS(v(ind_r3), delta_v_34, a_34(1), h, rhs_o2_34);
  v(ind_r3 + 1 : ind_r4) = A_34 \ C_34;  % Solving

  %=============================================================================
  % T r4 < r < r5
  %=============================================================================
  % NO
  % Neumann BC for zero mass flux
  u(end) = u(end - 1);
  R_NO_45 = para.lambda_t .* u(ind_r4 + 1 : end - 1);
  % RHS
  B_45 = r_coeff_no .* R_NO_45;
  B_45(1) = B_45(1) - (1 - a_45(1)) * u(ind_r4);
  B_45(end) = B_45(end) - (1 + a_45(end)) * u(end);
  u(ind_r4 + 1 : end - 1) = A_45 \ B_45;  % Solving
  % O2
  % Neumann BC for zero mass flux
  v(end) = v(end - 1);
  app_Km_45 = para.Km .* (1 + u(ind_r4 + 1 : end - 1) ./ gas_no.C_ref);
  R_O2_45 = gas_o2.Q_max_t .* v(ind_r4 + 1 : end - 1) ./...
      (v(ind_r4 + 1 : end - 1) + app_Km_45);
  % RHS
  C_45 = r_coeff_o2 .* R_O2_45;
  C_45(1) = C_45(1) - (1 - a_45(1)) * v(ind_r4);
  C_45(end) = C_45(end) - (1 + a_45(end)) * v(end);
  v(ind_r4 + 1 : end - 1) = A_45 \ C_45;  % Solving

  %=============================================================================
  % Steady state check here
  %=============================================================================
end

subplot(2, 1, 1);
plot(r(1 : ind_r1), u(1 : ind_r1),...
    r(ind_r1 + 1 : ind_r2), u(ind_r1 + 1 : ind_r2),...
    r(ind_r2 + 1 : ind_r3), u(ind_r2 + 1 : ind_r3),...
    r(ind_r3 + 1 : ind_r4), u(ind_r3 + 1 : ind_r4),...
    r(ind_r4 + 1 : end), u(ind_r4 + 1 : end));
legend('RBC', 'CFL', 'EC', 'VW', 'T');
subplot(2, 1, 2);
plot(r(1 : ind_r1), v(1 : ind_r1),...
    r(ind_r1 + 1 : ind_r2), v(ind_r1 + 1 : ind_r2),...
    r(ind_r2 + 1 : ind_r3), v(ind_r2 + 1 : ind_r3),...
    r(ind_r3 + 1 : ind_r4), v(ind_r3 + 1 : ind_r4),...
    r(ind_r4 + 1 : end), v(ind_r4 + 1 : end));
