%===============================================================================
% Creates an object with the general parameters of the model
% Parameters from the following paper:
% 1: Impact of the Fahraeus Effect on NO and O2 Biotransport A Computer Model
% 2: Two-dimensional transient model for prediction of arteriolar NO/O2
%    modulation by spatiotemporal variations in cell-free layer width
% \param offset_radius Offset radius of RBC core, in the current implementation
%        offset_radius is a percentage of normal_cfl
% \param offset_angle Offset angle of RBC core
% \param varargin Two cases: 1) Zero argument. 2) One argument.
%        1) When shape is a circle
%        2) \param deform Deformation of RBC core due to gravity effect
%===============================================================================
function obj = Parameters(num_coords, offset_radius, offset_angle, varargin)
  % General paramaters
  obj.alpha = 1.3;  % [uM/Torr], Solubility
  obj.int_r = 25.0;  % [um], Internal radius
  obj.Km = 1;  % [Torr], Michaelis constant in the absence of NO
  obj.lambda_b = 382.5;  % [1/s], Hb scavenging at 40% Hct
  obj.lambda_t = 1;  % [1/s], Tissue scavenging (from last slide)
  obj.lambda_vw = 1;  % [1/s], Vascular wall scavenging (from last slide)
  obj.len_EC = 2.5;  % [um], Endothelial cell width
  obj.len_T = 140;  % [um], Tissue layer width
  obj.len_VW = 10.0;  % [um], Vessel wall width (from paper)
  obj.R = obj.int_r + obj.len_EC + obj.len_T + obj.len_VW;
  % obj.wss = 1.5;  % [Pa], Wall shear stress
  wss_ref = 2.4;  % [Pa], Reference wall shear stress

  % CFL width calculation
  normal_cfl = (0.213 + 0.135) * obj.int_r / 2.0;
  rbc_core_radius = obj.int_r - normal_cfl;
  offset_radius = normal_cfl * offset_radius;
  switch length(varargin)
  case 0
    obj.shape = 'circle';
    quarter_coords = GetQuarterCoordinates(num_coords, offset_radius,...
        offset_angle, rbc_core_radius);
  case 1
    obj.shape = 'ellipse';
    obj.deform = varargin{1};
    quarter_coords = GetQuarterCoordinates(num_coords, offset_radius,...
        offset_angle, rbc_core_radius, rbc_core_radius, 0, obj.deform);
  otherwise
    error('Invalid shape');
  end
  obj.cfl = obj.int_r - quarter_coords;

  % Wall shear stress
  viscosity = 1.3;  % [cP], Plasma viscosity
  edge_vel = 1.0;  % [], Edge velocity
  wss = viscosity * edge_vel ./ obj.cfl;  % [Pa], Wall shear stress

  % Gas parameters
  obj.no.d_coeff = 3300.0;  % [um^2/s], Diffusion coefficient for NO
  obj.no.q_ref = 50.0;  % [uM/s], Reference/control NO production rate
  obj.no.C_ref = 27e-3;  % [uM], Reference NO concentration
  obj.no.R_max = wss / wss_ref * obj.no.q_ref;

  obj.o2.d_coeff = 2800.0;  % [um^2/s], Diffusion coefficient for O2
  obj.o2.Q_max_vw = 5.0;  % [uM/s], Max O2 consumption rate at vascular wall
  obj.o2.Q_max_t = 50.0;  % [uM/s], Max O2 consumption rate at tissue
  obj.o2.Km_eNOS = 4.7;  % [Torr]
  obj.o2.P = 70.0;  % [Torr], P_O2 in blood lumen
end
