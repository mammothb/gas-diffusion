%===============================================================================
% \param which_species 1 for NO, 2 for O2
%===============================================================================
function [EK, Ef] = MakeElemMat(para,...
    which_species,...
    compt_ind,...
    elem,...
    elem_node,...
    node_pos)
  % Set up Gauss points and weights for gauss quadrature
  gauss_pos = [0.5 - 1 / (2 * sqrt(3)), 0.5 + 1 / (2 * sqrt(3))];
  gauss_wgt = [0.5, 0.5];
  % Preallocate space for EK and Ef
  EK = zeros(para.nodes_per_elem, para.nodes_per_elem);
  Ef = zeros(para.nodes_per_elem, 1);
  ngp = length(gauss_pos);
  % Jacobian = abs(dx / dxi)
  J = abs(node_pos(elem_node(elem, 2)) - node_pos(elem_node(elem, 1)));
  for ii = 1 : ngp
    xi = gauss_pos(ii);
    % dPsi_n/dxi * dxi/dx = dPsi_n/dxi / J since it's 1D and the nodal position
    % scheme used
    B = [Psi(1, 1, xi), Psi(2, 1, xi)] ./ J;
    % Uses Gaussian quadrature to compute element stiffness matrix
    switch which_species
    case 1
      EK = EK + gauss_wgt(ii) * J .* (B' * B);
    case 2
      EK = EK + gauss_wgt(ii) * J .* (B' * B);
    otherwise
      error('Invalid species');
    end
  end
end
