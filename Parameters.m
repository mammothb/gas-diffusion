%===============================================================================
% Creates an object with the general parameters of the model
%===============================================================================
function obj = Parameters()
  obj.alpha = 1.3;  % [uM/Torr], Solubility
  obj.int_r = 25.0;  % [um], Internal radius
  obj.K_m = 1;  % [Torr], Michaelis constant in the absence of NO
  obj.K_m_app = 4.7;  % [Torr], Apparent Michaelis constant
  obj.lambda_b = 382.5;  % [1/s], Hb scavenging at 40% Hct
  obj.lambda_t = 0.1;  % [1/s], Tissue scavenging (from last slide)
  obj.lambda_vw = 0.1;  % [1/s], Vascular wall scavenging (from last slide)
  obj.len_EC = 2.5;  % [um], Endothelial cell width
  obj.len_T = 2500;  % [um], Tissue layer width
  obj.mu_p = 1.2;  % [cP], Plasma viscosity
  obj.n = 1.3;  % [1], Hill coefficient
  obj.R = obj.int_r + obj.len_EC + obj.len_T;
  obj.wss = 1.5;  % [Pa], Wall shear stress
  obj.wss_ref = 2.4;  % [Pa], Reference wall shear stress
end
