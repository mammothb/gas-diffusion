clear;
clf;

%% Model parameters
d_coeff = 1e-4;  % [cm^2/s], Diffusion coefficient
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
x = linspace(0, L, nx);
t = linspace(0, time, nt);
source = zeros(nx_i, 1);
source(q_start / h + 1 : q_end / h + 1) = q;
U = zeros(nx, nt);  % initial condition everywhere else is 0
U(1, :) = U_max * 0.01;  % boundary condition
U(end, :) = U_max * 0.01;  % boundary condition

delta = d_coeff * dt / h / h;  % Stability condition for explicit solution
I = speye(nx_i);
S = TriDiag(nx_i, -2, 1, 1);
A = I - delta / 2 * S;
B = I + delta / 2 * S;
C = zeros(nx_i, 1);

for tt = 1 : nt - 1
  for ii = 1 : nx_i
    C(ii) = 2 * dt * source(ii);  % add source term to RHS
  end
  % boundary condition
  C(1) = C(1) + delta / 2 * (U(1, tt) + U(1, tt + 1));
  C(end) = C(end) + delta / 2 * (U(end, tt) + U(end, tt + 1));
  % solving using crank-nicolson
  U(2 : end - 1, tt + 1) = A \ B * U(2 : end - 1, tt) + A \ C;
end

plot(x, U(:, 5 / dt + 1) / U_max);
hold on;
plot(x, U(:, 15 / dt + 1) / U_max);
hold on;
plot(x, U(:, 25 / dt + 1) / U_max);
hold on;
plot(x, U(:, 35 / dt + 1) / U_max);
% surf(t, x, U);
