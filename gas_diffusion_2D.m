clear;
clf;
%   y
%   ^
% H |
%   + - > x
%     L
U_exact = @(x,y,t) sin(pi * x) .* sin(pi * y) * exp(-2 * pi * pi * t);
% Model parameters
d_coeff = 1.0;  % [cm^2/s], Diffusion coefficient
L = 1.0;  % [cm], Length of channel
H = 1.0;  % [cm], Height of channel
T = 0.1;  % [s], Total simulation time
U_init = 0.1;  % [mM], Starting value for initial condition
% Dirichlet boundary values
U_N = 0;  %[mM], Concentration at north edge
U_S = 0;  %[mM], Concentration at south edge
U_E = 0;  %[mM], Concentration at east edge
U_W = 0;  %[mM], Concentration at west edge
% Simulation parameters
h = 0.01;  % [cm], space step
dt = h * h;  % [s], time step
nx = round(L / h) + 1;  % number of nodes in x direction
ny = round(H / h) + 1;  % number of nodes in y direction
nt = round(T / dt) + 1;  % number of time steps
nx_i = nx - 2;  % number of internal nodes in x direction
ny_i = ny - 2;  % number of internal nodes in y direction

% Initialization
x = linspace(0, L, nx);
y = linspace(0, H, ny);
t = linspace(0, T, nt);
U = zeros(ny, nx, nt);
% Initial condition
U_0 = zeros(ny, nx);
for jj = 1 : ny
  for ii = 1 : nx
    U_0(jj, ii) = sin(pi * ii * h) * sin(pi * jj * h);
  end
end
% Initial source term, for demo only
% U_0(floor(ny / 2) + 1, floor(nx / 2) + 1) = U_init * 2;
% Boundary conditions
% g0 = zeros(ny, 1);
% g1 = zeros(ny, 1);
% g2 = zeros(nx, 1);
% g3 = zeros(nx, 1);
g0 = ones(ny, 1) * U_W;
g1 = ones(ny, 1) * U_E;
g2 = ones(nx, 1) * U_S;
g3 = ones(nx, 1) * U_N;
% Boundary condition corner consistency
if g0(1) ~= g2(1) || g0(end) ~= g3(1) || g3(end) ~= g1(end) || g1(1) ~= g2(end)
  error('Boundary corner inconsistent');
end
for tt = 1 : nt
  U(:, 1, tt) = g0;
  U(:, end, tt) = g1;
  U(1, :, tt) = g2;
  U(end, :, tt) = g3;
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
U_vec = zeros(nx_i * ny_i, 1);
% U_new = U_old;
for jj = 1 : ny_i
  for ii = 1 : nx_i
    n = (jj - 1) * nx_i + ii;
    U_vec(n) = U_0(jj + 1, ii + 1);
  end
end
% Set up variables and matrix for C-N method
r = d_coeff * dt / h / h;  % stability condition for explicit solution
I = speye(nx_i);
A = TriDiag(nx_i, 2 * r, -r / 2, -r / 2);
Af = full(A);
b = zeros(nx_i, ny_i);
% combining compact forms
AA = BlockTriDiag(ny_i, A, -r / 2 * I, -r / 2 * I);
II = speye(nx_i * ny_i);
AAf = full(AA);
LL = (II + AA);
LLf = full(LL);
RR = (II - AA);
RRf = full(RR);
BB = zeros(nx_i * ny_i, 1);
% Start solving
for tt = 1 : nt - 1
  % maybe source term in here
  for jj = 1 : ny_i
    b(1, jj) = g0(jj + 1) + g0(jj + 1);
    b(end, jj) = g1(jj + 1) + g1(jj + 1);
  end
  for jj = 1 : ny_i
    ind1 = (jj - 1) * nx_i + 1;
    ind2 = jj * nx_i;
    % make necessary changes when boundary condition changes with time
    % extra terms for when slice is adjacent to boundary
    % maybe source term in here
    B2 = g2(2 : end - 1);
    B3 = g3(2 : end - 1);
    % Set up RHS
    BB(ind1 : ind2) = b(:, jj);
    if jj == 1
      BB(ind1 : ind2) = BB(ind1 : ind2) + B2 + B2;
    elseif jj == ny_i
      BB(ind1 : ind2) = BB(ind1 : ind2) + B3 + B3;
    end
  end
  U_vec = LL \ (RR * U_vec + BB);
  for jj = 1 : ny_i
    for ii = 1 : nx_i
      n = (jj - 1) * nx_i + ii;
      U(jj + 1, ii + 1, tt + 1) = U_vec(n);
    end
  end
end
% subplot(2, 2, 1);
% % contour(x, y, U(:, :, 1) / U_init);
% surf(x, y, U(:, :, 1) / U_init);
% subplot(2, 2, 2);
% % contour(x, y, U(:, :, round(0.02 / dt) + 1) / U_init);
% surf(x, y, U(:, :, round(0.02 / dt) + 1) / U_init);
% subplot(2, 2, 3);
% % contour(x, y, U(:, :, round(0.05 / dt) + 1) / U_init);
% surf(x, y, U(:, :, round(0.05 / dt) + 1) / U_init);
% subplot(2, 2, 4);
% % contour(x, y, U(:, :, round(0.07 / dt) + 1) / U_init);
% surf(x, y, U(:, :, round(0.07 / dt) + 1) / U_init);
disp(y(floor(ny / 2) + 1));
disp(t(round(0.07 / dt) + 1));
plot(x, U(floor(ny / 2) + 1, :, round(0.07 / dt) + 1));
hold on;
plot(x, U_exact(x, ones(1, ny) * H / 2, 0.07));
