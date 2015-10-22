clear;
clf;
%% Model parameters
para = Parameters();  % general parameters
gas_no = Gas('NO', para);  % create object for NO gas
gas_o2 = Gas('O2', para);  % create object for O2 gas
cfl = 3;  % [um], CFL width
lambda_core = para.lambda_b * 0.5 * (1 + para.int_r * para.int_r /...
    (para.int_r - cfl) / (para.int_r - cfl));

%% Simulation parameters
h = 0.5;  % [um], space step
if mod(para.R, h) > 1e-20
  error('Domain and space step incompatible');
end
nr = round(para.R / h) + 1;  % number of nodes in r direction
% nr_i = nr - 2;  % number of internal nodes
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
v = gas_o2.P * ones(nr, 1);  % solution for O2
r = linspace(0, para.R, nr);
r_01 = linspace(r(1), r(ind_r1), nr_01);
r_12 = linspace(r(ind_r1), r(ind_r2), nr_12);
r_23 = linspace(r(ind_r2), r(ind_r3), nr_23);
r_34 = linspace(r(ind_r3), r(end), nr_34);
r_45 = linspace(r(ind_r3), r(end), nr_45);
fprintf('r1: %.1f\nr2: %.1f\nr3: %.1f\n', r(ind_r1), r(ind_r2), r(ind_r3));

%% All R_species terms have been treated as if they are shifted to the RHS

for ii = 1 : 2
% RBC r0 < r < r1
% NO
a_01 = h / 2 ./ r(2 : ind_r1 - 1);
diags_01 = [[1 - a_01(2 : end)'; -2 * h; 0], [-2 * ones(nr_i_01 - 1, 1);...
    1.5 * h], [0; 1 + a_01(1 : end)']];
A_01 = spdiags(diags_01, [-1; 0; 1], nr_i_01, nr_i_01);
A_01(1, 1) = A_01(1, 1) - 1 + a_01(1);  % "Neumann" boundary at r = 0
A_01(end, end - 2) = h / 2;  % extra term for one-sided approximation of u'
% Neumann BC for zero mass flux
u(1) = u(2);
% Neumann BC for continuous mass flux
% Using one-side approximation
sigma_u_01 = (-3 * u(ind_r1) + 4 * u(ind_r1 + 1) - u(ind_r1 + 2)) / 2 / h;
% RHS
R_NO_01 = lambda_core .* u(2 : ind_r1 - 1);
B_01 = [h * h / gas_no.d_coeff .* R_NO_01; sigma_u_01];
% Solving
if ii ~= 1
  u(2:ind_r1) = A_01 \ B_01;
end
% O2
% Nothing to be done for O2 since P_o2 is constant in this region

%% CFL r1 < r < r2
% NO
a_12 = h / 2 ./ r(ind_r1 + 1 : ind_r2 - 1);
diags_12 = [[1 - a_12(2 : end)'; -2 * h; 0], [-2 * ones(nr_i_12 - 1, 1);...
    1.5 * h], [0; 1 + a_12(1 : end)']];
A_12 = spdiags(diags_12, [-1; 0; 1], nr_i_12, nr_i_12);
A_12(end, end - 2) = h / 2;
% Neumann BC for continuous mass flux
phi_u_12 = (-3 * u(ind_r2) + 4 * u(ind_r2 + 1) - u(ind_r2 + 2)) / 2 / h;
% RHS
B_12 = [-(1 - a_12(1)) * u(ind_r1); zeros(nr_i_12 - 2, 1); phi_u_12];
% Solving
if ii ~= 1
  u(ind_r1 + 1 : ind_r2) = A_12 \ B_12;
end
% O2
% Neumann BC for continous mass flux
phi_v_12 = (-3 * v(ind_r2) + 4 * v(ind_r2 + 1) - v(ind_r2 + 2)) / 2 / h;
% RHS
C_12 = [-(1 - a_12(1)) * v(ind_r1); zeros(nr_i_12 - 2, 1); phi_v_12];
% Solving
v(ind_r1 + 1 : ind_r2) = A_12 \ C_12;

%% EC r2 < r < r3
% NO
a_23 = h / 2 ./ r(ind_r2 + 1 : ind_r3 - 1);
% when using one-side approixmation for Neumann BC
diags_23 = [[1 - a_23(2 : end)'; -2 * h; 0], [-2 * ones(nr_i_23 - 1, 1);...
    1.5 * h], [0; 1 + a_23(1 : end)']];
A_23 = spdiags(diags_23, [-1; 0; 1], nr_i_23, nr_i_23);
A_23(end, end - 2) = h / 2;
% Neumann BC for continuous mass flux using one-sided approximation
gamma_u_23 = (-3 * u(ind_r3) + 4 * u(ind_r3 + 1) - u(ind_r3 + 2)) / 2 / h;
% NO production term
% one sided approx
R_NO_23 = -gas_no.R_max .* v(ind_r2 + 1 : ind_r3 - 1) ./...
    (v(ind_r2 + 1 : ind_r3 - 1) + gas_o2.Km_eNOS);
% RHS
B_23 = [h * h / gas_no.d_coeff .* R_NO_23; gamma_u_23];  % one sided approx
B_23(1) = B_23(1) - (1 - a_23(1)) * u(ind_r2);
% Solving
if ii ~= 1
  u(ind_r2 + 1 : ind_r3) = A_23 \ B_23;
end
% O2
% Neumann BC for continuous mass flux using one-sided approximation
gamma_v_23 = (-3 * v(ind_r3) + 4 * v(ind_r3 + 1) - v(ind_r3 + 2)) / 2 / h;
% RHS
% one sided approx
C_23 = [-R_NO_23 ./ gas_o2.d_coeff ./ para.alpha; gamma_v_23];
C_23(1) = C_23(1) - (1 - a_23(1)) * v(ind_r2);
% Solving
v(ind_r2 + 1 : ind_r3) = A_23 \ C_23;

%% VW r3 < r < r4
% NO
a_34 = h / 2 ./ r(ind_r3 + 1 : ind_r4 - 1);
diags_34 = [[1 - a_34(2 : end)'; -2 * h; 0], [-2 * ones(nr_i_34 - 1, 1);...
    1.5 * h], [0; 1 + a_34(1 : end)']];
A_34 = spdiags(diags_34, [-1; 0; 1], nr_i_34, nr_i_34);
A_34(end, end - 2) = h / 2;
% Neumann BC for continous mass flux
delta_u_34 = (-3 * u(ind_r4) + 4 * u(ind_r4 + 1) - u(ind_r4 + 2)) / 2 / h;
R_NO_34 = para.lambda_vw .* u(ind_r3 + 1 : ind_r4 - 1);
% RHS
B_34 = [h * h / gas_no.d_coeff .* R_NO_34; delta_u_34];
B_34(1) = B_34(1) - (1 - a_34(1)) * u(ind_r3);
% Solving
if ii ~= 1
  u(ind_r3 + 1 : ind_r4) = A_34 \ B_34;
end
% O2
% Neumann BC for continous mass flux
delta_v_34 = (-3 * v(ind_r4) + 4 * v(ind_r4 + 1) - v(ind_r4 + 2)) / 2 / h;
app_Km = para.Km * (1 + u(ind_r3 + 1 : ind_r4 - 1) ./ gas_no.C_ref);
R_O2_34 = gas_o2.Q_max_vw .* v(ind_r3 + 1 : ind_r4 - 1) ./...
    (v(ind_r3 + 1 : ind_r4 - 1) + app_Km);
% RHS
C_34 = [h * h / gas_o2.d_coeff / para.alpha .* R_O2_34; delta_v_34];
C_34(1) = C_34(1) - (1 - a_34(1)) * v(ind_r3);
% Solving
v(ind_r3 + 1 : ind_r4) = A_34 \ C_34;

%% T r4 < r < r5
% NO
a_45 = h / 2 ./ r(ind_r4 + 1 : end - 1);
diags_45 = [[1 - a_45(2 : end)'; 0], [-2 * ones(nr_i_45, 1)], [0;...
    1 + a_45(1 : end - 1)']];
A_45 = spdiags(diags_45, [-1; 0; 1], nr_i_45, nr_i_45);
A_45(end, end) = A_45(end, end) + 1 + a_45(end - 1);  % BC at r = r5
% Neumann BC for zero mass flux
u(end) = u(end - 1);
R_NO_45 = para.lambda_t .* u(ind_r4 + 1 : end - 1);
% RHS
B_45 = [h * h / gas_no.d_coeff .* R_NO_45];
B_45(1) = B_45(1) - (1 - a_45(1)) * u(ind_r4);
% Solving
if ii ~= 1
  u(ind_r4 + 1 : end - 1) = A_45 \ B_45;
end
% O2
% Neumann BC for zero mass flux
v(end) = v(end - 1);
app_Km = para.Km .* (1 + u(ind_r4 + 1 : end - 1) ./ gas_no.C_ref);
R_O2_45 = gas_o2.Q_max_t .* v(ind_r4 + 1 : end - 1) ./...
    (v(ind_r4 + 1 : end - 1) + app_Km);
% RHS
C_45 = [h * h / gas_o2.d_coeff / para.alpha .* R_O2_45];
C_45(1) = C_45(1) - (1 - a_45(1)) * v(ind_r4);
% Solving
v(ind_r4 + 1 : end - 1) = A_45 \ C_45;
end
subplot(2, 1, 1);
% plot(r(1 : ind_r4), u(1 : ind_r4));
plot(r, u);
subplot(2, 1, 2);
% plot(r(1 : ind_r4), v(1 : ind_r4));
plot(r, v);
