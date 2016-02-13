clf;
clear all;
H_d = 0.335;  % Discharge hematocrit
int_r = 27.1;  % [um], Internal radius
nr = 101;
r = linspace(0, int_r, nr);

% Guessing parameters
% Step 1: Initial guess for H_c
H_c = 0.9797 * H_d + 0.0404;  % Initial guess for core hematocrit
k_0 = 25.608 * H_c * H_c - 35.091 * H_c + 14.981;
k_inf = 1.3125 * H_c * H_c - 2.3001 * H_c + 2.6884;
gamma_c = 30.443 * H_c * H_c - 13.528 * H_c + 1.49;
a_0 = 1;
a_1 = 0;
a_2 = 0;

% Step 2: Calculate CFL width
cfl = a_2 * H_d * H_d + a_1 * H_d + a_0;
r1 = int_r - cfl;
shear_rate = zeros(1, nr);
J = 3732;
mu_p = 1.3;
% Step 3: Calculate shear rate
for ii = 1 : nr
  if r(ii) < r1
    core_sr = @(gamma) 0.5 * J / mu_p * r(ii) * (1.0 - 0.5 * H_c * (k_0 +...
        k_inf * sqrt(gamma / gamma_c)) / (1.0 + sqrt(gamma / gamma_c)))^2 -...
        gamma;
    shear_rate(ii) = fzero(core_sr, [0.0 100000.0]);
  else
    shear_rate(ii) = 0.5 * J * r(ii) / mu_p;
  end
end
% Step 4: Compute flow velocity
velocity = zeros(1, nr);
for ii = 1 : nr
  velocity(ii) = trapz(shear_rate(1 : ii));
end
subplot(2, 1, 1);
plot(r, shear_rate);
subplot(2, 1, 2);
plot(r, velocity);
