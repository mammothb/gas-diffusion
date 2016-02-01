function RevolvePlot(params, r0, theta0, varargin)
  switch length(varargin)
  case 0
    A = importdata('data_small.dat');
    radius = A(:, 1);
    data_width = (length(A(1, :)) - 1) / 2;
    u_start = 2;
    u_end = u_start + data_width - 1;
    v_start = u_end + 1;
    v_end = v_start + data_width - 1;
    u = A(:, u_start : u_end);
    v = A(:, v_start : v_end);
  case 3
    radius = varargin{1};
    u = varargin{2};
    v = varargin{3};
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
  len_coarse = length(u(1, :)) + 1;
  [theta_coarse, r_coarse] = meshgrid(linspace(0, 2 * pi, len_coarse), radius);
  [theta, r] = meshgrid(linspace(0, 2 * pi, nr), radius);
  [x, y] = pol2cart(theta, r);
  u_r_coarse = [u, u(:, 1)];
  u_r = interp2(theta_coarse, r_coarse, u_r_coarse, theta, r);
  v_r_coarse = [v, v(:, 1)];
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
  if strcmp(params.shape, 'circle')
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
  elseif strcmp(params.shape, 'ellipse')
    a = rbc_core_radius;
    b = rbc_core_radius;
    phi = 0;
    sqr_diff = b * b - a * a;
    sqr_sum = a * a + b * b;
    deform = params.deform;
    qtr_len = nr / 4;
    nr_1 = nr * 0.25 + 1;
    nr_2 = nr * 0.5;
    nr_3 = nr * 0.75 + 1;
    for ii = 1 : nr
      t = thetas(ii);
      R = sqr_diff * cos(2 * t - 2 * phi) + sqr_sum;
      P = r0 * (sqr_diff * cos(t + theta0 - 2 * phi) + sqr_sum *...
          cos(t - theta0));
      Q = sqrt(2) * a * b * sqrt(R - 2 * r0 * r0 * sin(t - theta0) *...
          sin(t - theta0));
      if ii <= nr_2
        deformation = 1.0 - deform * (1.0 - abs(ii - nr_1) / qtr_len);
      else
        deformation = 1.0 + deform * (1.0 - abs(ii - nr_3) / qtr_len);
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
  title('Predicted radial gradients for C_{NO}', 'FontName', 'Cambria Math');
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
  title('Predicted radial gradients for P_{O_2}', 'FontName', 'Cambria Math');
  xlabel('Distance away from the vessel center [\mum]', 'FontName',...
      'Cambria Math');
  ylabel('Distance away from the vessel center [\mum]', 'FontName',...
      'Cambria Math');
  zlabel('P_{O_2} [Torr]', 'FontName', 'Cambria Math');
end
