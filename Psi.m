%======================================================
% Function to calculate basis functions and derivatives
%
% num = local node number
% der = derivative: 0 = function; 1 = d/dxi_1
% xi = xi coordinate (range 0-1)
%======================================================
function result = Psi(num, der, xi)
  % Checks that input varialbes are valid
  if ~any(num == [1, 2]), error('Invalid local node number'); end
  if ~any(der == [0, 1]), error('Invalid derivative number'); end
  if xi > 1 || xi < 0, error('Invalid derivative number'); end

  values = [1 - xi, xi;...
                -1,  1];
  result = values(der + 1, num);
end
