clear;
close all;
% Read data. Data columns: r, u1, u2, u3, u4, u5, v1, v2, v3, v4, v5
% u is NO concentration, v is O2 concentration
A = importdata('data_small1.dat');
radius = A(:, 1);

% Find indexes of the interfaces
% r1_1 is r1 for CFL = 1
% r1_5 is r1 for CFL = 5
ind_r1_1 = find(24 == radius);
ind_r1_5 = find(20 == radius);
ind_r2 = find(25 == radius);
ind_r3 = find(27.5 == radius);
ind_r4 = find(37.5 == radius);

% Extract concentration values for plotting
u1 = A(:, 2);
u5 = A(:, 6);
v1 = A(:, 7);
v5 = A(:, 11);
nr = length(radius);
% Create polar coordinates
[theta, r] = meshgrid(linspace(0, 2 * pi, nr), radius);
[x, y] = pol2cart(theta, r);
u1_r = u1 * ones(1, nr);
u5_r = u5 * ones(1, nr);
v1_r = v1 * ones(1, nr);
v5_r = v5 * ones(1, nr);

% Initialize array for drawing lines at interface when CFL = 1
% Initialize at -1 so it will not overlap surface plot
layers_1 = -1 * ones(nr, nr);
layers_1(ind_r1_1, :) = ones(1, nr);
layers_1(ind_r2, :) = ones(1, nr);
layers_1(ind_r3, :) = ones(1, nr);
layers_1(ind_r4, :) = ones(1, nr);
% Set values at interface to be 71 so it is visible even in RBC core for O2
layers_1 = 71 .* layers_1;

% Interface lines when CFL = 5
layers_5 = -1 * ones(nr, nr);
layers_5(ind_r1_5, :) = ones(1, nr);
layers_5(ind_r2, :) = ones(1, nr);
layers_5(ind_r3, :) = ones(1, nr);
layers_5(ind_r4, :) = ones(1, nr);
layers_5 = 71 .* layers_5;

% u1
figure;
surf(x, y, u1_r, 'linestyle', 'none');
hold on;
surf(x, y, layers_1);
view(2);
colormap jet;
cb = colorbar;
ylabel(cb, 'C_{NO} [\muM]', 'FontName', 'Cambria Math', 'FontSize', 12);
caxis([0, 0.17]);
title('Predicted radial gradients for C_{NO} at \delta = 1', 'FontName', 'Cambria Math');
xlabel('Distance away from the vessel center [\mum]', 'FontName', 'Cambria Math');
ylabel('Distance away from the vessel center [\mum]', 'FontName', 'Cambria Math');
zlabel('C_{NO} [\muM]', 'FontName', 'Cambria Math');

% v1
figure;
surf(x, y, v1_r, 'linestyle', 'none');
hold on;
surf(x, y, layers_1);
view(2);
colormap jet;
cb = colorbar;
ylabel(cb, 'P_{O_2} [Torr]', 'FontName', 'Cambria Math', 'FontSize', 12);
caxis([0, 70]);
title('Predicted radial gradients for P_{O_2} at \delta = 1', 'FontName', 'Cambria Math');
xlabel('Distance away from the vessel center [\mum]', 'FontName', 'Cambria Math');
ylabel('Distance away from the vessel center [\mum]', 'FontName', 'Cambria Math');
zlabel('P_{O_2} [Torr]', 'FontName', 'Cambria Math');

% u5
figure;
surf(x, y, u5_r, 'linestyle', 'none');
hold on;
surf(x, y, layers_5);
view(2);
colormap jet;
cb = colorbar;
ylabel(cb, 'C_{NO} [\muM]', 'FontName', 'Cambria Math', 'FontSize', 12);
caxis([0, 0.17]);
title('Predicted radial gradients for C_{NO} at \delta = 5', 'FontName', 'Cambria Math');
xlabel('Distance away from the vessel center [\mum]', 'FontName', 'Cambria Math');
ylabel('Distance away from the vessel center [\mum]', 'FontName', 'Cambria Math');
zlabel('C_{NO} [\muM]', 'FontName', 'Cambria Math');

% v5
figure;
surf(x, y, v5_r, 'linestyle', 'none');
hold on;
surf(x, y, layers_5);
view(2);
colormap jet;
cb = colorbar;
ylabel(cb, 'P_{O_2} [Torr]', 'FontName', 'Cambria Math', 'FontSize', 12);
caxis([0, 70]);
title('Predicted radial gradients for P_{O_2} at \delta = 5', 'FontName', 'Cambria Math');
xlabel('Distance away from the vessel center [\mum]', 'FontName', 'Cambria Math');
ylabel('Distance away from the vessel center [\mum]', 'FontName', 'Cambria Math');
zlabel('P_{O_2} [Torr]', 'FontName', 'Cambria Math');