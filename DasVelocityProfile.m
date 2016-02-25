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

mesh(x, z, H);
xlabel('x');
ylabel('z');
zlabel('H');
zlim([0 0.8]);
