clf;
clear all;
H_d = 0.355;  % Discharge hematocrit
int_r = 27.1;  % [um], Internal radius
nr = 101; % Number of data
r = linspace(0, int_r, nr); % Create data from 0 to int_rad with 101 data
h = r(2) - r(1); % space-step
tolerance = 1e-6;
conv_flag = true;

% Guessing parameters
% Step 1: Initial guess for H_c
H_c = 0.9797 * H_d + 0.0404;  % Initial guess for core hematocrit (631)
k_0 = 25.608 * H_c * H_c - 35.091 * H_c + 14.981; % Figure 2 
k_inf = 1.3125 * H_c * H_c - 2.3001 * H_c + 2.6884; % Figure 2 
gamma_c = 30.443 * H_c * H_c - 13.528 * H_c + 1.49; % Figure 2 

% Step 2: Calculate CFL width
cfl = -2.265 * H_d * H_d - 1.4377 * H_d + 3.2131;  % Eqn 15
r1 = int_r - cfl; % Calculate RBC core
ind_r1 = round(r1 / h + 1);  % index for end of inner vessel

shear_rate = zeros(1, nr); % Create data of 1 row and 101 column
J = 3732; % Pressure Gradient, dyn/cm^3 (632)
mu_p = 1.3; % Plasma Viscosity 

% Step 3: Calculate shear rate
while conv_flag
  H_c_old = H_c;  % Copy old hematocrit for convergence test later
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
  
  % Step 5: Refine H_c
  H_c = H_d * trapz(velocity) / trapz(velocity(1 : ind_r1));
  
  % Step 6 Check that the new and old Hc is within the tolerance level
  abs_diff = abs(H_c - H_c_old);
  disp(abs_diff);
  if abs_diff < tolerance, conv_flag = false; end
end

subplot(2, 1, 1);
plot(r, shear_rate);
subplot(2, 1, 2);
plot(r, velocity);
