function RevolvePlot(params, r0, theta0, shape, varargin)
  switch length(varargin)
  case 0
    A = importdata('data_small.dat');
    radius = A(:, 1);
    u1 = A(:, 2);
    u2 = A(:, 3);
    u3 = A(:, 4);
    u4 = A(:, 5);
    v1 = A(:, 6);
    v2 = A(:, 7);
    v3 = A(:, 8);
    v4 = A(:, 9);
  case 3
    radius = varargin{1};
    u_ans = varargin{2};
    v_ans = varargin{3};
    u1 = u_ans(:, 1);
    u2 = u_ans(:, 2);
    u3 = u_ans(:, 3);
    u4 = u_ans(:, 4);
    v1 = v_ans(:, 1);
    v2 = v_ans(:, 2);
    v3 = v_ans(:, 3);
    v4 = v_ans(:, 4);
  end

  nr = length(radius);
  % Find indexes of the interfaces
  r2 = params.int_r;
  r3 = r2 + params.len_EC;
  r4 = r3 + params.len_VW;
  ind_r2 = find(r2 == radius);
  ind_r3 = find(r3 == radius);
  ind_r4 = find(r4 == radius);

  % Create polar coordinates
  [theta_coarse, r_coarse] = meshgrid(linspace(0, 2 * pi, 5), radius);
  [theta, r] = meshgrid(linspace(0, 2 * pi, nr), radius);
  [x, y] = pol2cart(theta, r);
  u_r_coarse = [u1, u2, u3, u4, u1];
  u_r = interp2(theta_coarse, r_coarse, u_r_coarse, theta, r);
  v_r_coarse = [v1, v2, v3, v4, v1];
  v_r = interp2(theta_coarse, r_coarse, v_r_coarse, theta, r);

  % Draw layers
  layers = -1 * ones(nr, nr);
  layers(ind_r2, :) = ones(1, nr);
  layers(ind_r3, :) = ones(1, nr);
  layers(ind_r4, :) = ones(1, nr);
  layers = 71 .* layers;

  % Drawing shifted RBC core
  thetas = linspace(0, 2 * pi, nr);
  normal_cfl = (0.213 + 0.135) * params.int_r / 2.0;
  rbc_core_radius = params.int_r - normal_cfl;
  r0 = normal_cfl * r0;
  if strcmp(shape, 'circle')
    r = rbc_core_radius;
    for ii = 1 : nr
      t = thetas(ii);
      r_front = r0 * cos(t - theta0);
      r_back = sqrt(r * r - r0 * r0 * sin(t - theta0) * sin(t - theta0));
      r_plus = r_front + r_back;
      r_minus = r_front - r_back;
      ind = find(round(r_plus, 1) == radius);
      layers(ind, ii) = 71;
    end
  elseif strcmp(shape, 'ellipse')
    a = rbc_core_radius;
    b = rbc_core_radius;
    phi = 0;
    sqr_diff = b * b - a * a;
    sqr_sum = a * a + b * b;
    deform = 0.1;
    nr_1 = nr * 0.25;
    nr_2 = nr * 0.5;
    nr_3 = nr * 0.75;
    for ii = 1 : nr
      t = thetas(ii);
      R = sqr_diff * cos(2 * t - 2 * phi) + sqr_sum;
      P = r0 * (sqr_diff * cos(t + theta0 - 2 * phi) + sqr_sum *...
          cos(t - theta0));
      Q = sqrt(2) * a * b * sqrt(R - 2 * r0 * r0 * sin(t - theta0) *...
          sin(t - theta0));
      if ii <= nr_2
        deformation = 1.0 + deform * (1.0 - abs(ii - nr_1) / nr_1);
      else
        deformation = 1.0 + deform * (1.0 - abs(ii - nr_3) / nr_1);
      end
      coord = (P + Q) / R * deformation;
      ind = find(round(coord, 1) == radius);
      layers(ind, ii) = 71;
    end
  end

  % u
  figure;
  surf(x, y, u_r, 'linestyle', 'none');
  hold on;
  surf(x, y, layers);
  view(2);
  colormap jet;
  cb = colorbar;
  title('Predicted radial gradients for C_{NO} at \delta = 1', 'FontName',...
      'Cambria Math');
  xlabel('Distance away from the vessel center [\mum]', 'FontName',...
      'Cambria Math');
  ylabel('Distance away from the vessel center [\mum]', 'FontName',...
      'Cambria Math');
  zlabel('C_{NO} [\muM]', 'FontName', 'Cambria Math');
  caxis([0, 0.07]);

  % v
  figure;
  surf(x, y, v_r, 'linestyle', 'none');
  hold on;
  surf(x, y, layers);
  view(2);
  colormap jet;
  cb = colorbar;
  ylabel(cb, 'P_{O_2} [Torr]', 'FontName', 'Cambria Math', 'FontSize', 12);
  caxis([0, 70]);
  title('Predicted radial gradients for P_{O_2} at \delta = 1', 'FontName',...
      'Cambria Math');
  xlabel('Distance away from the vessel center [\mum]', 'FontName',...
      'Cambria Math');
  ylabel('Distance away from the vessel center [\mum]', 'FontName',...
      'Cambria Math');
  zlabel('P_{O_2} [Torr]', 'FontName', 'Cambria Math');
end
