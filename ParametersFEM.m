function obj = ParametersFEM(cfl_width, varargin)
  % Default value
  h = 1;
  if nargin > 1, h = varargin{1}; end

  obj.alpha = 1.3;  % [uM/Torr], Solubility
  obj.int_r = 25.0;  % [um], Internal radius
  obj.Km = 1;  % [Torr], Michaelis constant in the absence of NO
  obj.lambda_b = 382.5;  % [1/s], Hb scavenging at 40% Hct
  obj.lambda_t = 1;  % [1/s], Tissue scavenging (from last slide)
  obj.lambda_vw = 1;  % [1/s], Vascular wall scavenging (from last slide)
  obj.len_EC = 2.5;  % [um], Endothelial cell width
  obj.len_T = 100;  % [um], Tissue layer width
  obj.len_VW = 10.0;  % [um], Vessel wall width (from paper)
  obj.R = obj.int_r + obj.len_EC + obj.len_T + obj.len_VW;
  obj.wss = 1.5;  % [Pa], Wall shear stress
  obj.wss_ref = 2.4;  % [Pa], Reference wall shear stress

  obj.no.d_coeff = 3300.0;  % [um^2/s], Diffusion coefficient for NO
  obj.no.q_ref = 50.0;  % [uM/s], Reference/control NO production rate
  obj.no.C_ref = 27e-3;  % [uM], Reference NO concentration
  obj.no.R_max = obj.wss / obj.wss_ref * obj.no.q_ref;

  obj.o2.d_coeff = 2800.0;  % [um^2/s], Diffusion coefficient for O2
  obj.o2.Q_max_vw = 5.0;  % [uM/s], Max O2 consumption rate at vascular wall
  obj.o2.Q_max_t = 50.0;  % [uM/s], Max O2 consumption rate at tissue
  obj.o2.Km_eNOS = 4.7;  % [Torr]
  obj.o2.P = 70.0;  % [Torr], P_O2 in blood lumen

  obj.lambda_core = obj.lambda_b / 2 * (1 + obj.int_r * obj.int_r /...
      (obj.int_r - cfl_width) * (obj.int_r - cfl_width));
  obj.r_coeff_no = h * h / obj.no.d_coeff;
  obj.r_coeff_o2 = h * h / obj.no.d_coeff / obj.alpha;

  obj.elem_per_compt = 10;  % Number of elements per compartment
  obj.num_compt = 5;  % Number of compartments
  obj.nodes_per_elem = 2;  % Number of (local) nodes per element

  obj.num_elems = obj.num_compt * obj.elem_per_compt;  % Total number of elements
  obj.num_nodes = obj.num_elems + 1;  % Total number of nodes

  % Different space steps in each compartment
  obj.h = zeros(obj.num_compt, 1);
  obj.h(1) = (obj.int_r - cfl_width) / obj.elem_per_compt;
  obj.h(2) = cfl_width / obj.elem_per_compt;
  obj.h(3) = obj.len_EC / obj.elem_per_compt;
  obj.h(4) = obj.len_VW / obj.elem_per_compt;
  obj.h(5) = obj.len_T / obj.elem_per_compt;
end
