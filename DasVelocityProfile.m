clf;
clear;
% Model parameters
% From paper
n = 20;  % Power law index
H_min = 0.0;  % Minimum hematocrit

% Arbitrary
H_max = 0.6;  % Maximum hematocrit
R1 = 25.0;  % Vessel radius
R2 = 20.0;  % RBC core radius
lambda = R2 / R1;  % Core thickness
h = 0.2 * R1;  % Minimum gap width
delta = h / R1;  % Dimensionless peripheral layer thickness
d = 0.1 * R1;  % Distance between two centers
epsilon = d / R1;  % Eccentricity

if epsilon > 1 - lambda, error('Eccentricity too high'); end

nn = 101;  % Number of nodes in each direction
grid_vector = linspace(-1, 1, nn);  % -1 and 1 are normalized values
[x, z] = meshgrid(grid_vector);
r = sqrt(x.^2 + z.^2);  % Radial coordinates
H = ones(nn, nn) .* H_min;  % Hematrocrit distribution

for ii = 1 : nn
  for jj = 1 : nn
    % if r(ii, jj) ~= 0
    %   rho = sqrt(epsilon^2 + r(ii, jj)^2 + 2 * r(ii, jj) * epsilon *...
    %       z(ii, jj) / r(ii, jj));
    % else
    %   rho = sqrt(epsilon^2 + r(ii, jj)^2);
    % end
    % %% Non-Axisymmetric hematocrit distribution
    % if rho < lambda
    %   rl = rho / lambda;  % For convenience
    %   H(ii, jj) = H_max - (H_max - H_min) * (0.5 * n * (n - 1) * rl^(n - 2) -...
    %       n * (n - 2) * rl^(n - 1) + 0.5 * (n - 1) * (n - 2) * rl^n);
    % end
    %% Axisymmetric hematocrit distribution
    if r(ii, jj) < 1 - delta
      H(ii, jj) = H_max - (H_max - H_min) * (0.5 * n * (n - 1) * (r(ii, jj) /...
      lambda)^(n - 2) - n * (n - 2) * (r(ii, jj) / lambda)^(n - 1) + 0.5 *...
      (n - 1) * (n - 2) * (r(ii, jj) / lambda) ^ n);
    end
  end
end

% Parameters for RK4
mu_p = 1.3;
P_g = 1.0;
I_e = 1.0;  % *
gamma_s = 1.0;  % *
% Model
k_0 = 0.275363 + 2 ./ (0.100158 + H);
k_inf = exp(1.3435 + H .* (-2.803 + H .* (2.711 - 0.6479 .* H)));
gamma_c = exp(-6.1508 + H .* (27.923 + H .* (-25.6 + 3.697 .* H)));
k = (k_0 + k_inf .* sqrt(I_e ./ gamma_c)) ./ (1 + sqrt(I_e ./ gamma_c));  % *
mu = mu_p ./ (1 - 0.5 .* k .* H) .^ 2;  % *
tau_w = 8 .* mu .* gamma_s;  % *
mu_inf = mu_p ./ (1 - 0.5 .* k_inf .* H) .^ 2;
Lambda = gamma_c .* ((1 - 0.5 .* k_0 .* H) ./ (1 - 0.5 .* k_inf .* H)) .^ 2;
tau_0 = mu_p .* gamma_c .* (0.5 .* H .* (k_0 - k_inf)) .^ 2 ./ (1 - 0.5 .*...
    k_inf .* H) .^ 4;
alpha = (sqrt(tau_0) + sqrt(mu_inf .* Lambda)) ./ sqrt(lambda .* tau_w);
q = (sqrt(tau_0) - sqrt(mu_inf .* Lambda)) ./ (sqrt(tau_0) + sqrt(mu_inf .*...
    Lambda));
C_r = sqrt(r - 2 .* alpha .* q .* sqrt(r) + alpha .^ 2);
% Arbitrary
h = 0.02;
R = 0 : h : 1;
incr = 1;
v = zeros(1, numel(R));
x_i = 50;
for ii = numel(R) - 1 : -1 : 1
  % Reverse sign when going from R1 to 0
  mu_inf_i = mu_inf(x_i, ii + x_i);
  r_i = r(x_i, ii + x_i);
  alpha_i = alpha(x_i, ii + x_i);ra
  q_i = q(x_i, ii + x_i);
  C_r_i = C_r(x_i, ii + x_i);
  dvdr = @(r_i, x) (P_g * R1 / 4 / mu_inf_i * (r_i - alpha_i * (1 + q_i) *...
      sqrt(r_i) + alpha_i^2 + (sqrt(r_i) - alpha_i) * C_r_i));
  k1 = h * dvdr(R(ii + incr), v(ii + incr));
  k2 = h * dvdr(R(ii + incr) + h / 2, v(ii + incr) + k1 / 2);
  k3 = h * dvdr(R(ii + incr) + h / 2, v(ii + incr) + k2 / 2);
  k4 = h * dvdr(R(ii + incr) + h, v(ii + incr) + k3);
  v(ii) = v(ii + incr) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
end
plot(R, v);

% mesh(x, z, H);
% xlabel('x');
% ylabel('z');
% zlabel('H');
% zlim([0 0.8]);
