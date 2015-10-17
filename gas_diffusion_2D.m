clear;
clf;
%   y
%   ^
% H |
%   + - > x
%     L
% Model parameters
d_coeff = 0.0001;  % [cm^2/s], Diffusion coefficient
L = 0.2;  % [cm], Length of channel
H = 0.2;  % [cm], Height of channel
T = 35;  % [s], Total simulation time
U_source = 1.0;
% Dirichlet boundary values
U_N = 0.0;  %[mM], Concentration at north edge
U_S = 0.0;  %[mM], Concentration at south edge
U_E = 0.0;  %[mM], Concentration at east edge
U_W = 0.0;  %[mM], Concentration at west edge
% Simulation parameters
h = 0.01;  % [cm], space step
dt = 0.05;  % [s], time step
nx = L / h + 1;  % number of nodes in x direction
ny = H / h + 1;  % number of nodes in y direction
nt = T / dt + 1;  % number of time steps
nx_i = nx - 2;  % number of internal nodes in x direction
ny_i = ny - 2;  % number of internal nodes in y direction
% Initialization
x = linspace(0, L, nx);
y = linspace(0, H, ny);
t = linspace(0, T, nt);
U = zeros(ny, nx, nt);
% Initial condition
U_0 = zeros(ny, nx);
% Initial source term, for demo only
U_0(floor(ny / 2) + 1, floor(nx / 2) + 1) = U_source;
% Boundary conditions
g0 = zeros(ny, 1);
g1 = zeros(ny, 1);
g2 = zeros(nx, 1);
g3 = zeros(nx, 1);
% Boundary condition corner consistency
if g0(1) ~= g2(1) || g0(end) ~= g3(1) || g3(end) ~= g1(end) || g1(1) ~= g2(end)
  error('Boundary corner inconsistent');
end
% Add BC to IC
for jj = 1 : ny
  U_0(jj, 1) = g0(jj);
  U_0(jj, end) = g1(jj);
end  % jj
for ii = 1 : nx
  U_0(1, ii) = g2(ii);
  U_0(end, ii) = g3(ii);
end  % ii
% Applying IC and BC
U(:, :, 1) = U_0;
% Set up temporary storage for solutions
U_old = zeros(nx_i * ny_i, 1);
U_new = U_old;
for jj = 1 : ny_i
  for ii = 1 : nx_i
    n = (jj - 1) * nx_i + ii;
    U_old(n) = U_0(jj + 1, ii + 1);
  end
end
% Set up variables and matrix for C-N method
r = d_coeff * dt / h / h;  % stability condition for explicit solution
I = speye(nx_i);
S = TriDiag(nx_i, 2 * r, -r / 2, -r / 2);
b = zeros(nx_i, 1);
% combinding compact forms
SS = BlockTriDiag(ny_i, S, -r / 2 * I, -r / 2 * I);
II = speye(nx_i * ny_i);
AA = (II + SS);
BB = (II - SS);
CC = zeros(nx_i * ny_i, 1);

for tt = 1 : nt - 1
  for jj = 1 : ny_i
    ind1 = (jj - 1) * nx_i + 1;
    ind2 = jj * nx_i;
    % make necessary changes when boundary condition changes with time
    % extra terms for when slice is adjacent to boundary
    B2 = g2(2 : end - 1);
    B3 = g3(2 : end - 1);
    % maybe source term in here
    b(1) = g0(jj) + g0(jj);
    b(end) = g1(jj) + g1(jj);
    % Set up RHS
    CC(ind1 : ind2) = b;
    if jj == 1
      CC(ind1 : ind2) = CC(ind1 : ind2) + B2 + B2;
    elseif jj == ny_i
      CC(ind1 : ind2) = CC(ind1 : ind2) + B3 + B3;
    end
  end
  U_new = AA \ BB * U_old + AA \ CC;
  for jj = 1 : ny_i
    for ii = 1 : nx_i
      n = (jj - 1) * nx_i + ii;
      U(jj + 1, ii + 1, tt + 1) = U_new(n);
    end
  end
  U_old = U_new;
end
subplot(2, 2, 1);
surf(x, y, U(:, :, 5 / dt + 1) / U_source);
% zlim([0, 0.1]);
subplot(2, 2, 2);
surf(x, y, U(:, :, 15 / dt + 1) / U_source);
% zlim([0, 0.1]);
subplot(2, 2, 3);
surf(x, y, U(:, :, 25 / dt + 1) / U_source);
% zlim([0, 0.1]);
subplot(2, 2, 4);
surf(x, y, U(:, :, 35 / dt + 1) / U_source);
% zlim([0, 0.1]);
