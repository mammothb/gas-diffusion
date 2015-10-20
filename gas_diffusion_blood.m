clear;
clf;
%% Model parameters
% General
alpha = 1.3;  % [uM/Torr], Solubility
lambda_b = 382.5;  % [1/s], Hb scavenging at 40% Hct
lambda_vw = 1.0;  % [1/s], Vascular wall scavenging (0.1?)
lambda_t = 1.0;  % [1/s], Tissue scavenging (0.1?)
K_m = 1;  % [Torr], Michaelis constant in the absence of NO
K_m_app = 4.7;  % [Torr], Apparent Michaelis constant
wss = 1.5;  % [Pa], Wall shear stress
mu_p = 1.2;  % [cP], Plasma viscosity
wss_ref = 2.4;  % [Pa], Reference wall shear stress
y_EC = 2.5;  % [um], Endothelial cell width
y_T = 2500;  % [um], Tissue layer width
n = 1.3;  % [1], Hill coefficient

% NO related stuff (using an object)
NO.D = 3300;  % [um^2/s], Diffusion coefficient for NO
NO.q_ref = 50;  % [uM/s], Reference/control NO production rate
NO.C_ref = 27e-3;  % [uM], Reference NO concentration

% Oxygen related stuff
O2.D = 2800;  % [um^2/s], Diffusion coefficient for O2
O2.Q_max_vw = 5;  % [uM/s], Max O2 consumption rate at vascular wall
O2.Q_max_t = 50;  % [uM/s], Max O2 consumption rate at tissue
O2.K_m_eNO = 4.7;  % [Torr]
O2.P = 70.0;  % [Torr], P_O2 in blood lumen

L = 0.2;  % [cm], Length of channel
time = 35;  % [s], Total simulation time
q = 0.5;  % [mM], Source term
q_start = 0.05;  % [cm], Source start point
q_len = 0.02;  % [cm], Length of source term
q_end = q_start + q_len;  % [cm], Source end point
% Boundary conditions
U_max = 20;  % [mM], Max concentration
% Simulation parameters
h = 0.001;  % [cm], space step
dt = 0.05;  % [s], time step
nx = L / h + 1;  % number of nodes in x direction
nx_i = nx - 2;  % number of internal nodes
nt = time / dt + 1;  % number of time steps
% Initialization
% x = linspace(0, L, nx);
% t = linspace(0, time, nt);
% source = zeros(nx_i, 1);
% source(q_start / h + 1 : q_end / h + 1) = q;
% U = zeros(nx, nt);  % initial condition everywhere else is 0
% U(1, :) = U_max * 0.01;  % boundary condition
% U(end, :) = U_max * 0.01;  % boundary condition
% % Set up variables and matrix for C-N method
% delta = d_coeff * dt / h / h;  % Stability condition for explicit solution
% I = speye(nx_i);
% S = TriDiag(nx_i, -2, 1, 1);
% A = I - delta / 2 * S;
% B = I + delta / 2 * S;
% C = zeros(nx_i, 1);
%
% for tt = 1 : nt - 1
%   for ii = 1 : nx_i
%     C(ii) = 2 * dt * source(ii);  % add source term to RHS
%   end
%   % boundary condition
%   C(1) = C(1) + delta / 2 * (U(1, tt) + U(1, tt + 1));
%   C(end) = C(end) + delta / 2 * (U(end, tt) + U(end, tt + 1));
%   % solving using crank-nicolson
%   U(2 : end - 1, tt + 1) = A \ B * U(2 : end - 1, tt) + A \ C;
% end
%
% plot(x, U(:, 5 / dt + 1) / U_max);
% hold on;
% plot(x, U(:, 15 / dt + 1) / U_max);
% hold on;
% plot(x, U(:, 25 / dt + 1) / U_max);
% hold on;
% plot(x, U(:, 35 / dt + 1) / U_max);
% % surf(t, x, U);
