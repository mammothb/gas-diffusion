%===============================================================================
% Creates an object with the general parameters of the model
% Parameters from the following paper:
% 1: Impact of the Fahraeus Effect on NO and O2 Biotransport A Computer Model
% 2: Two-dimensional transient model for prediction of arteriolar NO/O2
%    modulation by spatiotemporal variations in cell-free layer width
%===============================================================================
function obj = Parameters()
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
end
