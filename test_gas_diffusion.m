function tests = FunctionalityTests
  tests = functiontests(localfunctions);
end

function testMakeNOLHS_scheme2(testCase)
  scheme = 2;  % Select one-sided approximation scheme
  % Arbitrary values for variables to tests functionality
  omega = 1.5;  % Non-unity omega value
  nr = 10;  % Matrix dimensions
  a = linspace(2, nr - 1, nr - 2);  % Set vector 'a' such that value are 2 to 9
  loose_tol = 1e-6;  % Tolerance

  % Expected values, calculations in excel sheet
  % M_NO = 1 / omega * D + L
  % N_NO = (1 / omega - 1) * D - U
  expected_M_NO = [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                   -1, -1.333333333, 0, 0, 0, 0, 0, 0, 0, 0;
                   0, -2, -1.333333333, 0, 0, 0, 0, 0, 0, 0;
                   0, 0, -3, -1.333333333, 0, 0, 0, 0, 0, 0;
                   0, 0, 0, -4, -1.333333333, 0, 0, 0, 0, 0;
                   0, 0, 0, 0, -5, -1.333333333, 0, 0, 0, 0;
                   0, 0, 0, 0, 0, -6, -1.333333333, 0, 0, 0;
                   0, 0, 0, 0, 0, 0, -7, -1.333333333, 0, 0;
                   0, 0, 0, 0, 0, 0, 0, -8, -1.333333333, 0;
                   0, 0, 0, 0, 0, 0, 0, 1, -4, 2];
  expected_N_NO = [1, -4, 1, 0, 0, 0, 0, 0, 0, 0;
                   0, 0.666666667, -3, 0, 0, 0, 0, 0, 0, 0;
                   0, 0, 0.666666667, -4, 0, 0, 0, 0, 0, 0;
                   0, 0, 0, 0.666666667, -5, 0, 0, 0, 0, 0;
                   0, 0, 0, 0, 0.666666667, -6, 0, 0, 0, 0;
                   0, 0, 0, 0, 0, 0.666666667, -7, 0, 0, 0;
                   0, 0, 0, 0, 0, 0, 0.666666667, -8, 0, 0;
                   0, 0, 0, 0, 0, 0, 0, 0.666666667, -9, 0;
                   0, 0, 0, 0, 0, 0, 0, 0, 0.666666667, -10;
                   0, 0, 0, 0, 0, 0, 0, 0, 0, -1];

  % Convert to sparse to match output format of function
  expected_M_NO = sparse(expected_M_NO);
  expected_N_NO = sparse(expected_N_NO);

  % Output from function
  [M_NO, N_NO] = MakeNOLHS(scheme, omega, a, nr);

  % Test function output are close to expected output within the given tolerance
  verifyEqual(testCase, M_NO, expected_M_NO, 'RelTol', loose_tol);
  verifyEqual(testCase, N_NO, expected_N_NO, 'RelTol', loose_tol);
end

function testMakeO2LHS_scheme2(testCase)
  scheme = 2;  % Select one-sided approximation scheme
  % Arbitrary values for variables to tests functionality
  omega = 1.5;  % Non-unity omega value
  nr = 9;  % Matrix dimensions
  ind = 1;  % Index set to 1 to read the entirety of vector 'a'
  a = linspace(2, nr, nr - 1);  % Set vector 'a' such that value are 2 to 9
  loose_tol = 1e-6;  % Tolerance

  % Expected values, calculations in excel sheet
  % M_O2 = 1 / omega * D + L
  % N_O2 = (1 / omega - 1) * D - U
  expected_M_O2 = [0.666666667, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                   -1, -1.333333333, 0, 0, 0, 0, 0, 0, 0, 0;
                   0, -2, -1.333333333, 0, 0, 0, 0, 0, 0, 0;
                   0, 0, -3, -1.333333333, 0, 0, 0, 0, 0, 0;
                   0, 0, 0, -4, -1.333333333, 0, 0, 0, 0, 0;
                   0, 0, 0, 0, -5, -1.333333333, 0, 0, 0, 0;
                   0, 0, 0, 0, 0, -6, -1.333333333, 0, 0, 0;
                   0, 0, 0, 0, 0, 0, -7, -1.333333333, 0, 0;
                   0, 0, 0, 0, 0, 0, 0, -8, -1.333333333, 0;
                   0, 0, 0, 0, 0, 0, 0, 1, -4, 2];
  expected_N_O2 = [-0.333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                   0, 0.666666667, -3, 0, 0, 0, 0, 0, 0, 0;
                   0, 0, 0.666666667, -4, 0, 0, 0, 0, 0, 0;
                   0, 0, 0, 0.666666667, -5, 0, 0, 0, 0, 0;
                   0, 0, 0, 0, 0.666666667, -6, 0, 0, 0, 0;
                   0, 0, 0, 0, 0, 0.666666667, -7, 0, 0, 0;
                   0, 0, 0, 0, 0, 0, 0.666666667, -8, 0, 0;
                   0, 0, 0, 0, 0, 0, 0, 0.666666667, -9, 0;
                   0, 0, 0, 0, 0, 0, 0, 0, 0.666666667, -10;
                   0, 0, 0, 0, 0, 0, 0, 0, 0, -1];

  % Convert to sparse to match output format of function
  expected_M_O2 = sparse(expected_M_O2);
  expected_N_O2 = sparse(expected_N_O2);

  % Output from function
  [M_O2, N_O2] = MakeO2LHS(scheme, omega, a, nr, ind);

  % Test function output are close to expected output within the given tolerance
  verifyEqual(testCase, M_O2, expected_M_O2, 'RelTol', loose_tol);
  verifyEqual(testCase, N_O2, expected_N_O2, 'RelTol', loose_tol);
end

function testMakeNORHS_ConstantConcentration(testCase)
  offset_radius = 0;
  offset_angle = 0;
  shape = 'circle';
  index = 1;
  params = Parameters(offset_radius, offset_angle, shape);
  scheme = 2;
  loose_tol = 1e-6;  % Tolerance
  % Arbitrary values for variables to tests functionality
  r_coeff_no = 1.3;  % Coefficient for R_NO term
  lambda_core = 1.4;  % Scavenging rate in RBC core
  nr = 20;  % Total number of nodes
  r = linspace(1, nr, nr);  % r value at each nodes
  n_cmpt = 4;  % Number of nodes per compartment
  % Indexes to make the end of each compartment
  ind_r1 = n_cmpt;
  ind_r2 = n_cmpt * 2;
  ind_r3 = n_cmpt * 3;
  ind_r4 = n_cmpt * 4;
  % r values of each compartment
  r_01 = r(1 : ind_r1);
  r_23 = r(ind_r2 + 1 : ind_r3);
  r_34 = r(ind_r3 + 1 : ind_r4);
  r_45 = r(ind_r4 + 1 : end);
  % Number of nodes in CFL compartment
  nr_12 = n_cmpt;
  % NO and O2 values (Arbitrary)
  u = 1.1 * ones(nr, 1);
  v = 1.2 * ones(nr, 1);

  exp_F = [0; 2.002; 2.002; 2.002; 0; 0; 0; 0; -1.665352779; -1.665352779;...
  -1.665352779; -1.665352779; 1.43; 1.43; 1.43; 1.43; 1.43; 1.43; 1.43;...
      0];
  F = MakeNORHS(params, u, v, r_coeff_no, scheme, r_01, nr_12, r_23, r_34,...
      r_45, lambda_core, index);

  verifyEqual(testCase, F, exp_F, 'RelTol', loose_tol);
end

function testMakeNORHS_VaryingConcentration(testCase)
  offset_radius = 0;
  offset_angle = 0;
  shape = 'circle';
  index = 1;
  params = Parameters(offset_radius, offset_angle, shape);
  scheme = 2;
  loose_tol = 1e-6;  % Tolerance
  % Arbitrary values for variables to tests functionality
  r_coeff_no = 1.3;  % Coefficient for R_NO term
  lambda_core = 1.4;  % Scavenging rate in RBC core
  nr = 20;  % Total number of nodes
  r = linspace(1, nr, nr);  % r value at each nodes
  n_cmpt = 4;  % Number of nodes per compartment
  % Indexes to make the end of each compartment
  ind_r1 = n_cmpt;
  ind_r2 = n_cmpt * 2;
  ind_r3 = n_cmpt * 3;
  ind_r4 = n_cmpt * 4;
  % r values of each compartment
  r_01 = r(1 : ind_r1);
  r_23 = r(ind_r2 + 1 : ind_r3);
  r_34 = r(ind_r3 + 1 : ind_r4);
  r_45 = r(ind_r4 + 1 : end);
  % Number of nodes in CFL compartment
  nr_12 = n_cmpt;
  % NO and O2 values (Arbitrary)
  u = [1.1 : 0.1 : 3]';
  v = [3 : -0.1 : 1.1]';

  exp_F = [0; 2.184; 2.366; 2.548; 0; 0; 0; 0; -2.610661723; -2.528642271;...
      -2.444174476; -2.357147052; 2.99; 3.12; 3.25; 3.38; 3.51; 3.64; 3.77;...
      0];
  F = MakeNORHS(params, u, v, r_coeff_no, scheme, r_01, nr_12, r_23, r_34,...
      r_45, lambda_core, index);

  verifyEqual(testCase, F, exp_F, 'RelTol', loose_tol);
end

function testMakeO2RHS_ConstantConcentration(testCase)
  offset_radius = 0;
  offset_angle = 0;
  shape = 'circle';
  index = 1;
  params = Parameters(offset_radius, offset_angle, shape);
  scheme = 2;
  loose_tol = 1e-6;  % Tolerance
  % Arbitrary values for variables to tests functionality
  r_coeff_o2 = 1.3;  % Coefficient for R_O2 term
  nr = 20;  % Total number of nodes
  r = linspace(1, nr, nr);  % r value at each nodes
  n_cmpt = 4;  % Number of nodes per compartment
  % Indexes to make the end of each compartment
  ind_r2 = n_cmpt * 2;
  ind_r3 = n_cmpt * 3;
  ind_r4 = n_cmpt * 4;
  % r values of each compartment
  r_23 = r(ind_r2 + 1 : ind_r3);
  r_34 = r(ind_r3 + 1 : ind_r4);
  r_45 = r(ind_r4 + 1 : end);
  % Number of nodes in CFL compartment
  nr_12 = n_cmpt;
  % NO and O2 values (Arbitrary)
  u = 1.1 * ones(nr, 1);
  v = 1.2 * ones(nr, 1);

  exp_G = [70; 0; 0; 0; 0; 1.665352779; 1.665352779; 1.665352779;...
      1.665352779; 0.181645679; 0.181645679; 0.181645679; 0.181645679;...
      1.816456788; 1.816456788; 1.816456788; 0];
  G = MakeO2RHS(params, u, v, r_coeff_o2, scheme, nr_12, r_23, r_34, r_45,...
      index);

  verifyEqual(testCase, G, exp_G, 'RelTol', loose_tol);
end

function testMakeO2RHS_VaryingConcentration(testCase)
  offset_radius = 0;
  offset_angle = 0;
  shape = 'circle';
  index = 1;
  params = Parameters(offset_radius, offset_angle, shape);
  scheme = 2;
  loose_tol = 1e-6;  % Tolerance
  % Arbitrary values for variables to tests functionality
  r_coeff_o2 = 1.3;  % Coefficient for R_O2 term
  nr = 20;  % Total number of nodes
  r = linspace(1, nr, nr);  % r value at each nodes
  n_cmpt = 4;  % Number of nodes per compartment
  % Indexes to make the end of each compartment
  ind_r2 = n_cmpt * 2;
  ind_r3 = n_cmpt * 3;
  ind_r4 = n_cmpt * 4;
  % r values of each compartment
  r_23 = r(ind_r2 + 1 : ind_r3);
  r_34 = r(ind_r3 + 1 : ind_r4);
  r_45 = r(ind_r4 + 1 : end);
  % Number of nodes in CFL compartment
  nr_12 = n_cmpt;
  % NO and O2 values (Arbitrary)
  u = [1.1 : 0.1 : 3]';
  v = [3 : -0.1 : 1.1]';

  exp_G = [70; 0; 0; 0; 0; 2.610661723; 2.528642271; 2.444174476;...
      2.357147052; 0.132976932; 0.120647822; 0.109252198; 0.09868791;...
      0.888671875; 0.797141959; 0.711630736; 0];
  G = MakeO2RHS(params, u, v, r_coeff_o2, scheme, nr_12, r_23, r_34, r_45,...
      index);

  verifyEqual(testCase, G, exp_G, 'RelTol', loose_tol);
end

function testGetQuarterCoordinates_CircleZeroAngularOffset(testCase)
  radius = 10;  % Circle radius
  r0 = 5;  % Radial offset
  theta0 = 0;  % Angular offset
  loose_tol = 1e-6;  % Tolerance

  exp_coords = [15, 8.7, 5, 8.7];
  coords = GetQuarterCoordinates(r0, theta0, radius);

  verifyEqual(testCase, coords, exp_coords, 'RelTol', loose_tol);
end

function testGetQuarterCoordinates_CircleVariousAngularOffset(testCase)
  radius = 10;  % Circle radius
  r0 = 5;  % Radial offset
  theta0_1 = pi * 0.25;  % Angular offset 45 deg
  theta0_2 = pi * 0.75;  % Angular offset 135 deg
  theta0_3 = pi * 1.25;  % Angular offset 225 deg
  theta0_4 = pi * 1.75;  % Angular offset 315 deg
  theta0_5 = pi * 1.62;  % Arbitrary angular offset
  loose_tol = 1e-6;  % Tolerance

  exp_coords1 = [12.9, 12.9, 5.8, 5.8];
  exp_coords2 = [5.8, 12.9, 12.9, 5.8];
  exp_coords3 = [5.8, 5.8, 12.9, 12.9];
  exp_coords4 = [12.9, 5.8, 5.8, 12.9];
  exp_coords5 = [10.7, 5.2, 7, 14.5];

  coords1 = GetQuarterCoordinates(r0, theta0_1, radius);
  coords2 = GetQuarterCoordinates(r0, theta0_2, radius);
  coords3 = GetQuarterCoordinates(r0, theta0_3, radius);
  coords4 = GetQuarterCoordinates(r0, theta0_4, radius);
  coords5 = GetQuarterCoordinates(r0, theta0_5, radius);

  verifyEqual(testCase, coords1, exp_coords1, 'RelTol', loose_tol);
  verifyEqual(testCase, coords2, exp_coords2, 'RelTol', loose_tol);
  verifyEqual(testCase, coords3, exp_coords3, 'RelTol', loose_tol);
  verifyEqual(testCase, coords4, exp_coords4, 'RelTol', loose_tol);
  verifyEqual(testCase, coords5, exp_coords5, 'RelTol', loose_tol);
end

function testGetQuarterCoordinates_EllipseZeroAngularOffset(testCase)
  a = 8;  % Ellipse major axis length
  b = 6;  % Ellipse minor axis length
  r0 = 5;  % Radial offset
  theta0 = 0;  % Angular offset
  phi = 0;  % Ellipse rotation
  loose_tol = 1e-6;  % Tolerance

  exp_coords = [13, 4.7, 3, 4.7];
  coords = GetQuarterCoordinates(r0, theta0, a, b, phi);

  verifyEqual(testCase, coords, exp_coords, 'RelTol', loose_tol);
end

function testGetQuarterCoordinates_EllipseVariousAngularOffset(testCase)
  a = 8;  % Ellipse major axis length
  b = 6;  % Ellipse minor axis length
  r0 = 5;  % Radial offset
  theta0_1 = pi * 0.25;  % Angular offset 45 deg
  theta0_2 = pi * 0.75;  % Angular offset 135 deg
  theta0_3 = pi * 1.25;  % Angular offset 225 deg
  theta0_4 = pi * 1.75;  % Angular offset 315 deg
  theta0_5 = pi * 1.62;  % Arbitrary angular offset
  phi = 0;  % Ellipse rotation
  loose_tol = 1e-6;  % Tolerance

  exp_coords1 = [10, 8.9, 2.9, 1.8];
  exp_coords2 = [2.9, 8.9, 10, 1.8];
  exp_coords3 = [2.9, 1.8, 10, 8.9];
  exp_coords4 = [10, 1.8, 2.9, 8.9];
  exp_coords5 = [6.9, 1.2, 3.2, 10.5];

  coords1 = GetQuarterCoordinates(r0, theta0_1, a, b, phi);
  coords2 = GetQuarterCoordinates(r0, theta0_2, a, b, phi);
  coords3 = GetQuarterCoordinates(r0, theta0_3, a, b, phi);
  coords4 = GetQuarterCoordinates(r0, theta0_4, a, b, phi);
  coords5 = GetQuarterCoordinates(r0, theta0_5, a, b, phi);

  verifyEqual(testCase, coords1, exp_coords1, 'RelTol', loose_tol);
  verifyEqual(testCase, coords2, exp_coords2, 'RelTol', loose_tol);
  verifyEqual(testCase, coords3, exp_coords3, 'RelTol', loose_tol);
  verifyEqual(testCase, coords4, exp_coords4, 'RelTol', loose_tol);
  verifyEqual(testCase, coords5, exp_coords5, 'RelTol', loose_tol);
end

function testGetQuarterCoordinates_EllipseVariousRotation(testCase)
  a = 8;  % Ellipse major axis length
  b = 6;  % Ellipse minor axis length
  r0 = 5;  % Radial offset
  theta0 = pi * 1.62;  % Arbitrary angular offset
  phi_1 = pi * 0.25;  % Rotation 45 deg
  phi_2 = pi * 0.75;  % Rotation 135 deg
  phi_3 = pi * 1.25;  % Rotation 225 deg
  phi_4 = pi * 1.75;  % Rotation 315 deg
  phi_5 = pi * 1.62;  % Arbitrary rotation

  loose_tol = 1e-6;  % Tolerance

  exp_coords1 = [8.3, 1.4, 2, 11.7];
  exp_coords2 = [5.7, 2.4, 4.6, 10.7];
  exp_coords3 = [8.3, 1.4, 2, 11.7];
  exp_coords4 = [5.7, 2.4, 4.6, 10.7];
  exp_coords5 = [6.1, 3.1, 3.9, 11.5];

  coords1 = GetQuarterCoordinates(r0, theta0, a, b, phi_1);
  coords2 = GetQuarterCoordinates(r0, theta0, a, b, phi_2);
  coords3 = GetQuarterCoordinates(r0, theta0, a, b, phi_3);
  coords4 = GetQuarterCoordinates(r0, theta0, a, b, phi_4);
  coords5 = GetQuarterCoordinates(r0, theta0, a, b, phi_5);

  verifyEqual(testCase, coords1, exp_coords1, 'RelTol', loose_tol);
  verifyEqual(testCase, coords2, exp_coords2, 'RelTol', loose_tol);
  verifyEqual(testCase, coords3, exp_coords3, 'RelTol', loose_tol);
  verifyEqual(testCase, coords4, exp_coords4, 'RelTol', loose_tol);
  verifyEqual(testCase, coords5, exp_coords5, 'RelTol', loose_tol);
end
