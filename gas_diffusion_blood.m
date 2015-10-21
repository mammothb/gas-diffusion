clear;
clf;
%% Model parameters
para = Parameters();  % general parameters
gas_no = Gas('NO');  % create object for NO gas
gas_o2 = Gas('O2');  % create object for O2 gas
cfl = 5.0;  % [um], CFL width
lambda_core = para.lambda_b * 0.5 * (1 + para.int_r * para.int_r /...
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
nr_01 = ind_r1;  % number of nodes for r0 < r < r1
nr_12 = ind_r2 - ind_r1;  % number of nodes for r1 < r < r2
nr_23 = ind_r3 - ind_r2;  % number of nodes for r2 < r < r3
nr_35 = nr - ind_r3;  % number of nodes for r3 < r < r5
nr_i_01 = nr_01 - 1;  % size of matrix for r0 < r < r1
nr_i_12 = nr_12;  % size of matrix for r1 < r < r2
nr_i_23 = nr_23;  % size of matrix for r2 < r < r3
nr_i_35 = nr_35 - 1;  % size of matrix for r3 < r < r5

%% Initialization
u = zeros(nr, 1);  % solution for NO
v = zeros(nr, 1);  % solution for O2
r = linspace(0, para.R, nr);
r_01 = linspace(r(1), r(ind_r1), nr_01);
r_12 = linspace(r(ind_r1 + 1), r(ind_r2), nr_12);
r_23 = linspace(r(ind_r2 + 1), r(ind_r3), nr_23);
r_35 = linspace(r(ind_r3 + 1), r(end), nr_35);
fprintf('r1: %.1f\nr2: %.1f\nr3: %.1f\n', r(ind_r1), r(ind_r2), r(ind_r3));

%% RBC
a_01 = h / 2 ./ r(2 : ind_r1 - 1);
diags_01 = [[1 - a_01(2 : end)'; -2 * h; 0], [-2 * ones(nr_i_01 - 1, 1);...
    1.5 * h], [0; 1 + a_01(1 : end)']];
A_01 = spdiags(diags_01, [-1; 0; 1], nr_i_01, nr_i_01);
A_01(1, 1) = A_01(1, 1) - (1 - a_01(1));  % Neumann boundary at r = 0
A_01(end, end - 2) = h / 2;  % extra term for one-sided approximation of u'
% Neumann BC for continuous mass flux
sigma_01 = (-3 * u(ind_r1) + 4 * u(ind_r1 + 1) - u(ind_r1 + 2)) / 2 / h;
% RHS
B_01 = [-lambda_core * u(2 : ind_r1 - 1); sigma_01];
% Solving RBC portion
u(2:ind_r1) = A_01 \ B_01;

%% CFL
a_12 = h / 2 ./ r(ind_r1 + 2 : ind_r2 - 1);
aa = [0; 1 - a_12(2 : end)'; -2 * h; 0];
bb = [-1.5 * h; -2 * ones(nr_i_12 - 2, 1); 1.5 * h];
cc = [0; 2 * h; 1 + a_12(1 : end)'];
diags_12 = [[0; 1 - a_12(2 : end)'; -2 * h; 0], [-1.5 * h;...
    -2 * ones(nr_i_12 - 2, 1); 1.5 * h], [0; 2 * h; 1 + a_12(1 : end)']];
A_12 = spdiags(diags_12, [-1; 0; 1], nr_i_12, nr_i_12);
A_12(1, 3) = -h / 2;
A_12(end, end - 2) = h / 2;
% Neumann BC for continuous mass flux
sigma_12 = (3 * u(ind_r1) - 4 * u(ind_r1 - 1) + u(ind_r1 - 2)) / 2 / h;
phi_12 = (-3 * u(ind_r2) + 4 * u(ind_r2 + 1) - u(ind_r2 + 2)) / 2 / h;
B_12 = [sigma_12; zeros(nr_i_12 - 2, 1); phi_12];
u(ind_r1 + 1 : ind_r2) = A_12 \ B_12;
plot(r, u);

% P = h / 2 ./ r(2 : end - 1);
% diags = [[1 - P(2 : end)'; 0], -2 * ones(nr_i, 1), [0; 1 + P(1 : end - 1)']];
% A = spdiags(diags, [-1; 0; 1], nr_i, nr_i);
% Af = full(A);

% Neumann boundary conditions
u(1) = u(2);  % Zero mass flux at r = 0
u(end) = u(end - 1);  % Zero mass flux at r = R
% main = spalloc(nr_i, nr_i, 3 * nr_i);
