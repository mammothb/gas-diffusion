%===============================================================================
% Creates an object with the general parameters of the model
% Parameters from the following paper:
% 1: Impact of the Fahraeus Effect on NO and O2 Biotransport A Computer Model
% 2: Two-dimensional transient model for prediction of arteriolar NO/O2
%    modulation by spatiotemporal variations in cell-free layer width
%===============================================================================
function obj = Parameters()
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
  obj.wss = 1.5;  % [Pa], Wall shear stress
  obj.wss_ref = 2.4;  % [Pa], Reference wall shear stress

  % Gas parameters
  obj.no.d_coeff = 3300.0;  % [um^2/s], Diffusion coefficient for NO
  obj.no.q_ref = 50.0;  % [uM/s], Reference/control NO production rate
  obj.no.C_ref = 27e-3;  % [uM], Reference NO concentration
  obj.no.R_max = obj.wss / obj.wss_ref * obj.no.q_ref;

  obj.o2.d_coeff = 2800.0;  % [um^2/s], Diffusion coefficient for O2
  obj.o2.Q_max_vw = 5.0;  % [uM/s], Max O2 consumption rate at vascular wall
  obj.o2.Q_max_t = 50.0;  % [uM/s], Max O2 consumption rate at tissue
  obj.o2.Km_eNOS = 4.7;  % [Torr]
  obj.o2.P = 70.0;  % [Torr], P_O2 in blood lumen
end
