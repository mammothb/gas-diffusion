%===============================================================================
% Creates an object for various gas types NO and O2
% \param which_gas indicates which gas to create
% \param param parameters object
%===============================================================================
function obj = Gas(which_gas, param)
  if which_gas == 'NO'
    obj.d_coeff = 3300.0;  % [um^2/s], Diffusion coefficient for NO
    obj.q_ref = 50.0;  % [uM/s], Reference/control NO production rate
    obj.C_ref = 27e-3;  % [uM], Reference NO concentration
    obj.R_max = param.wss / param.wss_ref * obj.q_ref;
  elseif which_gas == 'O2'
    obj.d_coeff = 2800.0;  % [um^2/s], Diffusion coefficient for O2
    obj.Q_max_vw = 5.0;  % [uM/s], Max O2 consumption rate at vascular wall
    obj.Q_max_t = 50.0;  % [uM/s], Max O2 consumption rate at tissue
    obj.Km_eNOS = 4.7;  % [Torr]
    obj.P = 70.0;  % [Torr], P_O2 in blood lumen
  end
end
