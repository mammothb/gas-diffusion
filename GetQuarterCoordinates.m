%===============================================================================
% This function calculates the coordinates at the four points where the off-
% center RBC core intersects with main axes of the blood vessel lumen
%
% \param radius Radius of the RBC core
% \param r0 The radial offset of the RBC core
% \param phi The angular offset of the RBC core (in radians)
% \return Coordinates at the four points, labeled in the following manner
%           - 2 -
%         /       \
%        3         1
%         \       /
%           - 4 -
%===============================================================================
function coordinates = GetQuarterCoordinates(radius, r0, phi)
  x_offset = r0 * cos(phi);
  y_offset = r0 * sin(phi);
  x_dist = sqrt(radius * radius + x_offset * x_offset);
  y_dist = sqrt(radius * radius + y_offset * y_offset);

  coordinates(1) = y_dist + x_offset;
  coordinates(2) = x_dist + y_offset;
  coordinates(3) = y_dist - x_offset;
  coordinates(4) = x_dist - y_offset;

  % Round up coordinates to nearest 1dp
  coordinates = round(coordinates, 1);
end
