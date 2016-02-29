clf;
clear;
% Model parameters
% From paper
n = 20;  % Power law index
H_min = 0.0;  % Minimum hematocrit
H_t = 0.25;  % Tube hematorcrit
% Arbitrary
R1 = 25.0;  % Vessel radius
R2 = 0.7 * R1;  % RBC core radius
lambda = R2 / R1;  % Core thickness
h = 0.2 * R1;  % Minimum gap width
delta = h / R1;  % Dimensionless peripheral layer thickness
d = 0.0 * R1;  % Distance between two centers
epsilon = d / R1;  % Eccentricity

% Only axisymmetric
H_max = H_min + (H_t - H_min) * (n + 1) * (n + 2) / (n - 1) / (n - 2) /...
    lambda ^ 2;

if epsilon > 1 - lambda, error('Eccentricity too high'); end

nn = 101;  % Number of nodes in each direction
grid_vector = linspace(-1, 1, nn);  % -1 and 1 are normalized values
[x, z] = meshgrid(grid_vector);
r = sqrt(x.^2 + z.^2);  % Radial coordinates
H = ones(nn, nn) .* H_min;  % Hematrocrit distribution

for ii = 1 : nn
  for jj = 1 : nn
    if r(ii, jj) ~= 0
      rho = sqrt(epsilon^2 + r(ii, jj)^2 + 2 * r(ii, jj) * epsilon *...
          z(ii, jj) / r(ii, jj));
    else
      rho = sqrt(epsilon^2 + r(ii, jj)^2);
    end
    %% Non-Axisymmetric hematocrit distribution
    if rho < lambda
      rl = rho / lambda;  % For convenience
      H(ii, jj) = H_max - (H_max - H_min) * (0.5 * n * (n - 1) * rl^(n - 2) -...
          n * (n - 2) * rl^(n - 1) + 0.5 * (n - 1) * (n - 2) * rl^n);
    end
    %% Axisymmetric hematocrit distribution
    % if r(ii, jj) < 1 - delta
    %   H(ii, jj) = H_max - (H_max - H_min) * (0.5 * n * (n - 1) * (r(ii, jj) /...
    %   lambda)^(n - 2) - n * (n - 2) * (r(ii, jj) / lambda)^(n - 1) + 0.5 *...
    %   (n - 1) * (n - 2) * (r(ii, jj) / lambda) ^ n);
    % end
  end
end

% Parameters for RK4
fmin = 0;
fmax = 5000;
tau_wm = 1.4;  % dyne / cm^2, Characteristic wall shear stress
mu_p = 1.3;
% P_g = 2 * tau_wm * 1e-5 * 1e8 / R1;
P_g = 4;
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
mu_inf_m = max(max(mu_inf));
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
% for ii = numel(R) - 1 : -1 : 1
%   x_i = 51;
%   z_i = x_i + ii;
%   % Reverse sign when going from R1 to 0
%   mu_inf_i = mu_inf(z_i, x_i);
%   r_i = r(z_i, x_i);
%   alpha_i = alpha(z_i, x_i);
%   q_i = q(z_i, x_i);
%   C_r_i = C_r(z_i, x_i);
%   dvdr = @(r_i, x) P_g * R1 / 4 / mu_inf_i * (r_i - alpha_i * (1 + q_i) *...
%       sqrt(r_i) + alpha_i^2 + (sqrt(r_i) - alpha_i) * C_r_i);
%   k1 = h * dvdr(R(ii + incr), v(ii + incr));
%   k2 = h * dvdr(R(ii + incr) + h / 2, v(ii + incr) + k1 / 2);
%   k3 = h * dvdr(R(ii + incr) + h / 2, v(ii + incr) + k2 / 2);
%   k4 = h * dvdr(R(ii + incr) + h, v(ii + incr) + k3);
%   v(ii) = v(ii + incr) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
% end
for ii = numel(R) - 1 : -1 : 1
  x_i = 51;
  z_i = x_i + ii;
  % Model Parameters
  k_0_i = k_0(z_i, x_i);
  k_inf_i = k_inf(z_i, x_i);
  gamma_c_i = gamma_c(z_i, x_i);
  H_i = H(z_i, x_i);
  tau_0_i = tau_0(z_i, x_i);
  mu_inf_i = mu_inf(z_i, x_i);
  Lambda_i = Lambda(z_i, x_i);
  q_i = q(z_i, x_i);
  r_i = r(z_i, x_i);
  k = @(I_e) (k_0_i + k_inf_i * sqrt(I_e / gamma_c_i)) / (1 + sqrt(I_e /...
      gamma_c_i));  % *
  mu = @(I_e) mu_p / (1 - 0.5 * k(I_e) * H_i) ^ 2 / mu_inf_m;  % *
  tau_w = @(I_e) 8 * mu(I_e) * gamma_s;  % *
  alpha_i = @(I_e) (sqrt(tau_0_i) + sqrt(mu_inf_i * Lambda_i)) /...
      sqrt(lambda * tau_w(I_e));
  C_r_i = @(I_e) sqrt(r_i - 2 * alpha_i(I_e) * q_i * sqrt(r_i) +...
      alpha_i(I_e) ^ 2);
  % Reverse sign when going from R1 to 0
  r_i = r(z_i, x_i);
  dvdr = @(r_i, x) P_g * R1 / 4 / mu_inf_i * (r_i - alpha_i * (1 + q_i) *...
      sqrt(r_i) + alpha_i^2 + (sqrt(r_i) - alpha_i) * C_r_i);
  r1 = R(ii + incr);
  r2 = R(ii + incr) + h / 2;
  r3 = R(ii + incr) + h / 2;
  r4 = R(ii + incr) + h;
  dvdr1 = @(I_e) P_g * R1 / 4 / mu_inf_i * (r1 - alpha_i(I_e) * (1 + q_i) *...
      sqrt(r1) + alpha_i(I_e)^2 + (sqrt(r1) - alpha_i(I_e)) * C_r_i(I_e)) - I_e;
  dvdr2 = @(I_e) P_g * R1 / 4 / mu_inf_i * (r2 - alpha_i(I_e) * (1 + q_i) *...
      sqrt(r2) + alpha_i(I_e)^2 + (sqrt(r2) - alpha_i(I_e)) * C_r_i(I_e)) - I_e;
  dvdr3 = @(I_e) P_g * R1 / 4 / mu_inf_i * (r3 - alpha_i(I_e) * (1 + q_i) *...
      sqrt(r3) + alpha_i(I_e)^2 + (sqrt(r3) - alpha_i(I_e)) * C_r_i(I_e)) - I_e;
  dvdr4 = @(I_e) P_g * R1 / 4 / mu_inf_i * (r4 - alpha_i(I_e) * (1 + q_i) *...
      sqrt(r4) + alpha_i(I_e)^2 + (sqrt(r4) - alpha_i(I_e)) * C_r_i(I_e)) - I_e;
  k1 = h * fzero(dvdr1, [fmin fmax]);
  k2 = h * fzero(dvdr2, [fmin fmax]);
  k3 = h * fzero(dvdr3, [fmin fmax]);
  k4 = h * fzero(dvdr4, [fmin fmax]);
  % k1 = h * dvdr(R(ii + incr), v(ii + incr));
  % k2 = h * dvdr(R(ii + incr) + h / 2, v(ii + incr) + k1 / 2);
  % k3 = h * dvdr(R(ii + incr) + h / 2, v(ii + incr) + k2 / 2);
  % k4 = h * dvdr(R(ii + incr) + h, v(ii + incr) + k3);
  v(ii) = v(ii + incr) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
end

plot(R, v * tau_wm * R1 / 2 / mu_inf_m);

% mesh(x, z, H);
% xlabel('x');
% ylabel('z');
% zlabel('H');
% zlim([0 0.8]);
