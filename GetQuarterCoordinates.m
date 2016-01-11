%===============================================================================
% This function calculates the coordinates at the four points where the off-
% center RBC core intersects with main axes of the blood vessel lumen
%
% \param r0 The radial offset of the RBC core
% \param theta0 The angular offset of the RBC core (in radians)
% \param varargin Two cases: 1) One argument. 2) Three arguments.
%        1) \param radius Radius of RBC core
%        2) \param a Length along major axis
%           \param b Length along minor axis
%           \param phi Rotation of RBC core along its center
% \return Coordinates at the four points, labeled in the following manner
%           - 2 -
%         /       \
%        3         1
%         \       /
%           - 4 -
%===============================================================================
function coordinates = GetQuarterCoordinates(r0, theta0, varargin)
  thetas = [0, 0.5 * pi, pi, 1.5 * pi];
  num_coords = length(thetas);
  coordinates = zeros(1, num_coords);

  switch length(varargin)
  case 1
    radius = varargin{1};
    for ii = 1 : num_coords
      theta = thetas(ii);
      r_front = r0 * cos(theta - theta0);
      r_back = sqrt(radius * radius - r0 * r0 * sin(theta - theta0) *...
          sin(theta - theta0));
      % r_front + r_back sometimes only describes half the circle, use
      % r_front - r_back when that occurs, have not occurred yet based on tests
      % performed
      r_plus = r_front + r_back;
      r_minus = r_front - r_back;
      if any(abs(coordinates - r_plus) < 0.001 * r0) && r_minus > 0
        coordinates(ii) = r_minus;
      else
        coordinates(ii) = r_plus;
      end
    end
  case 3
    a = varargin{1};
    b = varargin{2};
    phi = varargin{3};
    sqr_diff = b * b - a * a;
    sqr_sum = a * a + b * b;
    for ii = 1 : num_coords
      theta = thetas(ii);
      R = sqr_diff * cos(2 * theta - 2 * phi) + sqr_sum;
      P = r0 * (sqr_diff * cos(theta + theta0 - 2 * phi) + sqr_sum *...
          cos(theta - theta0));
      disp(P);
      Q = sqrt(2) * a * b * sqrt(R - 2 * r0 * r0 * sin(theta - theta0) *...
          sin(theta - theta0));
      coordinates(ii) = (P + Q) / R;
    end
  otherwise
    error('Incorrect number of input arguments');
  end

  % Round up coordinates to nearest 1dp
  coordinates = round(coordinates, 1);
end
