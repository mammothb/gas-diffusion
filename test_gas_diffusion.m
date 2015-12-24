function tests = FunctionalityTests
  tests = functiontests(localfunctions);
end

function testMakeNOLHS_scheme2(testCase)
  % Arbitrary values for variables to tests functionality
  scheme = 2;  % Select one-sided approximation scheme
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
  % Arbitrary values for variables to tests functionality
  scheme = 2;  % Select one-sided approximation scheme
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

function testGetQuarterCoordinates_ZeroAngularOffset(testCase)
  radius = 10;  % Circle radius
  r0 = 5;  % Radial offset
  phi = 0;  % Angular offset
  loose_tol = 1e-6;  % Tolerance

  exp_coords = [15, 11.2, 5, 11.2];
  coords = GetQuarterCoordinates(radius, r0, phi);

  verifyEqual(testCase, coords, exp_coords, 'RelTol', loose_tol);
end

function testGetQuarterCoordinates_VariousAngularOffset(testCase)
  radius = 10;  % Circle radius
  r0 = 5;  % Radial offset
  phi1 = pi * 0.25;  % Angular offset 45 deg
  phi2 = pi * 0.75;  % Angular offset 135 deg
  phi3 = pi * 1.25;  % Angular offset 225 deg
  phi4 = pi * 1.75;  % Angular offset 315 deg
  phi5 = pi * 1.62;  % Arbitrary angular offset
  loose_tol = 1e-6;  % Tolerance

  exp_coords1 = [14.1, 14.1, 7.1, 7.1];
  exp_coords2 = [7.1, 14.1, 14.1, 7.1];
  exp_coords3 = [7.1, 7.1, 14.1, 14.1];
  exp_coords4 = [14.1, 7.1, 7.1, 14.1];
  exp_coords5 = [12.9, 5.5, 9.2, 14.8];

  coords1 = GetQuarterCoordinates(radius, r0, phi1);
  coords2 = GetQuarterCoordinates(radius, r0, phi2);
  coords3 = GetQuarterCoordinates(radius, r0, phi3);
  coords4 = GetQuarterCoordinates(radius, r0, phi4);
  coords5 = GetQuarterCoordinates(radius, r0, phi5);

  verifyEqual(testCase, coords1, exp_coords1, 'RelTol', loose_tol);
  verifyEqual(testCase, coords2, exp_coords2, 'RelTol', loose_tol);
  verifyEqual(testCase, coords3, exp_coords3, 'RelTol', loose_tol);
  verifyEqual(testCase, coords4, exp_coords4, 'RelTol', loose_tol);
  verifyEqual(testCase, coords5, exp_coords5, 'RelTol', loose_tol);
end
