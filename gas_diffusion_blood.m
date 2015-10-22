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
nr_i_01 = nr_01 - 1;  % size of matrix for r0 < r < r1
nr_i_12 = nr_12;  % size of matrix for r1 < r < r2
nr_i_23 = nr_23;  % size of matrix for r2 < r < r3
nr_i_34 = nr_34;  % size of matrix for r3 < r < r4
nr_i_45 = nr_45 - 1;  % size of matrix for r3 < r < r4

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

% RBC r0 < r < r1
% NO
% a_01 = h / 2 ./ r(2 : ind_r1 - 1);

for ii = 1 : 2
%% LHS
% NO
a = h / 2 ./ r(2 : end - 1);
diags = [[1 - a(2 : end)'; 0], -2 * ones(nr_i, 1), [0; 1 + a(1 : end - 1)']];
A = spdiags(diags, [-1; 0; 1], nr_i, nr_i);
A(1, 1) = A(1, 1) - 1 + a(1);  % "Neumann" boundary at r = 0
A(end, end) = A(end, end) + 1 + a(end);  % "Neumann" boundary at r = 0
% O2
a_o2 = h / 2 ./ r(ind_r1 + 1 : end - 1);
% aaa = [zeros(nr_i_01, 1); 1 - a_o2(2 : end)'; 0];
% bbb = [ones(nr_i_01, 1); -2 * ones(nr_i - nr_i_01, 1)];
% ccc = [0; zeros(nr_i_01, 1); 1 + a_o2(1 : end - 1)'];
diags_o2 = [[zeros(nr_i_01, 1); 1 - a_o2(2 : end)'; 0],...
    [ones(nr_i_01, 1); -2 * ones(nr_i - nr_i_01, 1)],...
    [0; zeros(nr_i_01, 1); 1 + a_o2(1 : end - 1)']];
A_o2 = spdiags(diags_o2, [-1; 0; 1], nr_i, nr_i);

%%% RHS
B = zeros(nr_i, 1);
C = zeros(nr_i, 1);
%% RBC r0 < r < r1
% NO
B_01 = h * h * lambda_core / gas_no.d_coeff * u(2 : ind_r1);
% O2
% set to 70?
C_01 = gas_o2.P * ones(nr_i_01, 1);

%% CFL r1 < r < r2
% do nothing since both are 0
B_12 = zeros(nr_12 - 1, 1);
C_12 = zeros(nr_12 - 1, 1);

%% EC r2 < r < r3
% NO
R_NO_23 = h * h * gas_no.R_max .* v(ind_r2 + 1 : ind_r3) ./ (v(ind_r2 + 1 :...
    ind_r3) + gas_o2.K_m_eNOS);
B_23 = -R_NO_23 ./ gas_no.d_coeff;
% O2
C_23 = -R_NO_23 ./ gas_o2.d_coeff ./ para.alpha;

%% VW r3 < r < r4
% NO
B_34 = h * h / gas_no.d_coeff * para.lambda_vw .* u(ind_r3 + 1 : ind_r4);
% O2
app_Km = para.K_m .* (1 + u(ind_r3 + 1 : ind_r4) ./ gas_no.C_ref);
C_34 = h * h / gas_o2.d_coeff / para.alpha * gas_o2.Q_max_vw .*...
    v(ind_r3 + 1 : ind_r4) ./ (v(ind_r3 + 1 : ind_r4) + app_Km);

%% T r4 < r < r5
% NO
B_45 = h * h / gas_no.d_coeff * para.lambda_t .* u(ind_r4 + 1 : end - 1);
% O2
app_Km = para.K_m .* (1 + u(ind_r4 + 1 : end - 1) ./ gas_no.C_ref);
C_45 = h * h / gas_o2.d_coeff / para.alpha * gas_o2.Q_max_vw .*...
    v(ind_r4 + 1 : end - 1) ./ (v(ind_r4 + 1 : end - 1) + app_Km);

BB = [B_01; B_12; B_23; B_34; B_45];
u(2 : end - 1) = A \ BB;
CC = [C_01; C_12; C_23; C_34; C_45];
v(2 : end - 1) = A_o2 \ CC;
end

% Af = full(A);
% aa = det(Af);
% disp(aa);
% diags_01 = [[1 - a_01(2 : end)'; -2 * h; 0], [-2 * ones(nr_i_01 - 1, 1);...
%     1.5 * h], [0; 1 + a_01(1 : end)']];
% A_01 = spdiags(diags_01, [-1; 0; 1], nr_i_01, nr_i_01);
% A_01(1, 1) = A_01(1, 1) - 1 + a_01(1);  % Neumann boundary at r = 0
% A_01(end, end - 2) = h / 2;  % extra term for one-sided approximation of u'
% % Neumann BC for zero mass flux
% u(1) = u(2);
% % Neumann BC for continuous mass flux
% % Using one-side approximation
% sigma_u_01 = (-3 * u(ind_r1) + 4 * u(ind_r1 + 1) - u(ind_r1 + 2)) / 2 / h;
% % % Using central difference and extra point
% % sigma_u_01 = (u(ind_r1 + 1) - u(ind_r1 - 1)) / 2 / h;
% % RHS
% B_01 = [-lambda_core / gas_no.d_coeff * u(2 : ind_r1 - 1); sigma_u_01];
% % Solving
% % u(2:ind_r1) = A_01 \ B_01;
% % O2
% % Nothing to be done for O2 since P_o2 is constant in this region

% %% CFL r1 < r < r2
% % NO
% a_12 = h / 2 ./ r(ind_r1 + 1 : ind_r2 - 1);
% diags_12 = [[1 - a_12(1 : end)'; -2 * h; 0], [-1.5 * h;...
%     -2 * ones(nr_i_12 - 2, 1); 1.5 * h], [0; 2 * h; 1 + a_12(1 : end)']];
% A_12 = spdiags(diags_12, [-1; 0; 1], nr_i_12, nr_i_12);
% A_12(1, 3) = -h / 2;
% A_12(end, end - 2) = h / 2;
% % Neumann BC for continuous mass flux
% sigma_u_12 = (3 * u(ind_r1) - 4 * u(ind_r1 - 1) + u(ind_r1 - 2)) / 2 / h;
% phi_u_12 = (-3 * u(ind_r2) + 4 * u(ind_r2 + 1) - u(ind_r2 + 2)) / 2 / h;
% % RHS
% B_12 = [sigma_u_12; zeros(nr_i_12 - 2, 1); phi_u_12];
% % Solving
% % u(ind_r1 : ind_r2) = A_12 \ B_12;
% % O2
% % Neumann BC for continous mass flux
% sigma_v_12 = (3 * v(ind_r1) - 4 * v(ind_r1 - 1) + v(ind_r1 - 2)) / 2 / h;
% phi_v_12 = (-3 * v(ind_r2) + 4 * v(ind_r2 + 1) - v(ind_r2 + 2)) / 2 / h;
% % RHS
% C_12 = [sigma_v_12; zeros(nr_i_12 - 2, 1); phi_v_12];
% % Solving
% v(ind_r1 : ind_r2) = A_12 \ C_12;

%% EC r2 < r < r3
% NO
% a_23 = h / 2 ./ r(ind_r2 + 1 : ind_r3 - 1);
% % when using one-side approixmation for Neumann BC
% diags_23 = [[1 - a_23(1 : end)'; -2 * h; 0], [-1.5 * h;...
%     -2 * ones(nr_i_23 - 2, 1); 1.5 * h], [0; 2 * h; 1 + a_23(1 : end)']];
% A_23 = spdiags(diags_23, [-1; 0; 1], nr_i_23, nr_i_23);
% A_23(1, 3) = -h / 2;
% A_23(end, end - 2) = h / 2;
% when using extra node and centra difference for Neuman BC
% diags_23 = [[1 - a_23(1 : end)'; 1; 0], [-1;...
%     -2 * ones(nr_i_23 - 2, 1); -1], [0; 1; 1 + a_23(1 : end)']];
% A_23 = spdiags(diags_23, [-1; 0; 1], nr_i_23, nr_i_23);
% Neumann BC for continuous mass flux using one-sided approximation
% phi_u_23 = (3 * u(ind_r2) - 4 * u(ind_r2 - 1) + u(ind_r2 - 2)) / 2 / h;
% gamma_u_23 = (-3 * u(ind_r3) + 4 * u(ind_r3 + 1) - u(ind_r3 + 2)) / 2 / h;
% Neumann BC for continuous mass flux using extra node and central difference
% phi_u_23 = (u(ind_r2 + 1) - u(ind_r2 - 1)) / 2 / h;
% gamma_u_23 = (u(ind_r3 + 1) - u(ind_r3 - 1)) / 2 / h;
% NO production term
% one sided approx
% R_NO_23 = h * h * gas_no.R_max / gas_no.d_coeff .* v(ind_r2 + 1 :...
%     ind_r3 - 1) ./ (v(ind_r2 + 1 : ind_r3 - 1) + gas_o2.K_m_eNOS);
% extra node and central approx
% R_NO_23 = h * h * gas_no.R_max / gas_no.d_coeff .* v(ind_r2 : ind_r3) ./...
%     (v(ind_r2 : ind_r3) + gas_o2.K_m_eNOS);
% RHS
% B_23 = [phi_u_23; -R_NO_23; gamma_u_23];  % one sided approx
% B_23 = -R_NO_23;  % extra node and central approx

% Solving
% u(ind_r2 : ind_r3) = A_23 \ B_23;
% O2
% % Neumann BC for continuous mass flux using one-sided approximation
% phi_v_23 = (3 * v(ind_r2) - 4 * v(ind_r2 - 1) + v(ind_r2 - 2)) / 2 / h;
% gamma_v_23 = (-3 * v(ind_r3) + 4 * v(ind_r3 + 1) - v(ind_r3 + 2)) / 2 / h;
% Neumann BC for continuous mass flux using extra node and central difference
% phi_v_23 = (v(ind_r2 + 1) - v(ind_r2 - 1)) / 2 / h;
% gamma_v_23 = (v(ind_r3 + 1) - v(ind_r3 - 1)) / 2 / h;
% R_O2_23 = R_NO_23 .* (gas_no.d_coeff / gas_o2.d_coeff / para.alpha);
% RHS
% C_23 = [phi_v_23; -R_O2_23; gamma_v_23];  % one sided approx
% C_23 = [phi_v_23; -R_O2_23; gamma_v_23];  % extra node and central approx
% Solving
% v(ind_r2 : ind_r3) = A_23 \ C_23;

% %% VW r3 < r < r4
% % NO
% a_34 = h / 2 ./ r(ind_r3 + 1 : ind_r4 - 1);
% diags_34 = [[1 - a_34(1 : end)'; -2 * h; 0], [-1.5 * h;...
%     -2 * ones(nr_i_34 - 2, 1); 1.5 * h], [0; 2 * h; 1 + a_34(1 : end)']];
% A_34 = spdiags(diags_34, [-1; 0; 1], nr_i_34, nr_i_34);
% A_34(1, 3) = -h / 2;
% A_34(end, end - 2) = h / 2;
% % Neumann BC for continous mass flux
% gamma_u_34 = (3 * u(ind_r3) - 4 * u(ind_r3 - 1) + u(ind_r3 - 2)) / 2 / h;
% delta_u_34 = (-3 * u(ind_r4) + 4 * u(ind_r4 + 1) - u(ind_r4 + 2)) / 2 / h;
% R_NO_34 = h * h / gas_no.d_coeff * para.lambda_t * u(ind_r3 + 1 : ind_r4 - 1);
% % RHS
% B_34 = [gamma_u_34; R_NO_34; delta_u_34];
% % Solving
% % u(ind_r3 : ind_r4) = A_34 \ B_34;
% % O2
% % Neumann BC for continous mass flux
% gamma_v_34 = (3 * v(ind_r3) - 4 * v(ind_r3 - 1) + v(ind_r3 - 2)) / 2 / h;
% delta_v_34 = (-3 * v(ind_r4) + 4 * v(ind_r4 + 1) - v(ind_r4 + 2)) / 2 / h;
% app_Km = para.K_m * (1 + u(ind_r3 + 1 : ind_r4 - 1) ./ gas_no.C_ref);
% R_O2_34 = h * h / gas_o2.d_coeff * gas_o2.Q_max_vw .* v(ind_r3 + 1 :...
%     ind_r4 - 1) ./ (v(ind_r3 + 1 : ind_r4 - 1) + app_Km);
% % RHS
% C_34 = [gamma_v_34; R_O2_34; delta_v_34];
% % Solving
% v(ind_r3 : ind_r4) = A_34 \ C_34;

% %% T r4 < r < r5
% % NO
% a_45 = h / 2 ./ r(ind_r4 + 1 : end - 1);
% diags_45 = [[1 - a_45(1 : end)'; 0], [-1.5 * h; -2 * ones(nr_i_45 - 1, 1)],...
%     [0; 2 * h; 1 + a_45(1 : end - 1)']];
% A_45 = spdiags(diags_45, [-1; 0; 1], nr_i_45, nr_i_45);
% A_45(1, 3) = -h / 2;
% A_45(end, end) = A_45(end, end) + 1 + a_45(end - 1);  % BC at r = r5
% % Neumann BC for continous mass flux
% delta_u_45 = (3 * u(ind_r4) - 4 * u(ind_r4 - 1) + u(ind_r4 - 2)) / 2 / h;
% % Neumann BC for zero mass flux
% u(end) = u(end - 1);
% R_NO_45 = h * h / gas_no.d_coeff * para.lambda_t * u(ind_r4 + 1 : end - 1);
% % RHS
% B_45 = [delta_u_45; R_NO_45];
% % Solving
% % u(ind_r4 : end - 1) = A_45 \ B_45;
% % O2
% % Neumann BC for continous mass flux
% delta_v_45 = (3 * v(ind_r4) - 4 * v(ind_r4 - 1) + v(ind_r4 - 2)) / 2 / h;
% % Neumann BC for zero mass flux
% v(end) = v(end - 1);
% app_Km = para.K_m * (1 + u(ind_r4 + 1 : end - 1) ./ gas_no.C_ref);
% R_O2_45 = h * h / gas_o2.d_coeff * gas_o2.Q_max_t .* v(ind_r4 + 1 :...
%     end - 1) ./ (v(ind_r4 + 1 : end - 1) + app_Km);
% % RHS
% C_45 = [delta_v_45; R_O2_45];
% % Solving
% v(ind_r4 : end - 1) = A_45 \ C_45;

subplot(2, 1, 1);
% plot(r(ind_r2 + 1 : ind_r3), u(ind_r2 + 1 : ind_r3));
plot(r, u);
subplot(2, 1, 2);
plot(r, v);
