{ Fill in the following sections (removing comment marks ! if necessary),
  and delete those that are unused.}
TITLE 'Gas diffusion'     { the problem identification }
COORDINATES cylinder1  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  u              { NO concentration }
  v              { O2 partial pressure }
SELECT         { method controls }
  ERRLIM = 1e-6
DEFINITIONS    { parameter definitions }
  D_u = 3300	 { NO diffusion coefficient }
  D_v = 2800	 { O2 diffusion coefficient }
  alpha = 1.3	 { O2 solubility }
  R_u = 0		 { NO reaction term }
  R_v = 0        { O2 reaction term }
  int_r = 25	 { Inner radius }
  cfl = 5		 { CFL width }
  ec = 2.5       { EC width }
  vw = 10        { VW width }
  t = 100        { T width }
  Km_eNOS = 4.7  { K_{m,eNOS} for O2 }
  q_NO_ref = 50	 { Reference NO production }
  tau = 1.5
  tau_ref = 2.4
  R_NO_max = q_NO_ref * tau / tau_ref	 { Max NO production rate }
  Q_O2_max_vw = 5
  Q_O2_max_t = 50
  lambda_b = 382.5
  lambda_core = lambda_b * ((int_r - cfl)^2 + int_r^2) / (2 * (int_r - cfl)^2)
  lambda_vw = 1
  lambda_t = 1
  Km = 1
  C_ref = 27e-3
  r1 = int_r - cfl
  r2 = int_r
  r3 = r2 + ec
  r4 = r3 + vw
  r5 = r4 + t
! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  u: div(grad(u)) = R_u / D_u { one possibility }
  v: div(grad(v)) = R_v / D_v / alpha { one possibility }
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { RBC }
	R_u = lambda_core * u
	R_v = 0
    START(0)
	point natural(u) = 0  { Zero flux at the center }
	point natural(v) = 0  { Zero flux at the center }
    LINE TO (r1)
	point value(v) = 70
  REGION 2       { CFL }
	R_u = 0
	R_v = 0
    START(r1)
    LINE TO (r2)
  REGION 3       { EC }
	R_u = -R_NO_max * v / (v + Km_eNOS)
	R_v = -R_u
    START(r2)
    LINE TO (r3)
  REGION 4       { VW }
	R_u = lambda_vw * u
    R_v = Q_O2_max_vw * v / (v + Km * (1 + u / C_ref))
    START(r3)
    LINE TO (r4)
  REGION 5       { T }
	R_u = lambda_t * u
	R_v = Q_O2_max_t * v / (v + Km * (1 + u / C_ref))
    START(r4)
    LINE TO (r5)
	point natural(u) = 0
	point natural(v) = 0
MONITORS         { show progress }
PLOTS            { save result displays }
  ELEVATION(u) FROM (0) to (r5)
  ELEVATION(v) FROM (0) to (r5)
  TABLE(R, u, v) export format "#2"
END
