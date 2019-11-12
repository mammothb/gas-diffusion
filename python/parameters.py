class Parameters:
    # Physical parameters
    def __init__(self):
        self.alpha = 1.3  # [uM/Torr], Solubility
        self.lambda_b = 382.5  # [1/s], Hb scavenging at 40# Hct
        self.lambda_t = 1  # [1/s], Tissue scavenging (from last slide)
        self.lambda_vw = 1  # [1/s], Vascular wall scavenging (from last slide)
        self.Km = 1  # [Torr], Michaelis constant in the absence of NO

        self.int_r = 25  # [um], Internal radius
        self.len_EC = 2.5  # [um], Endothelial cell width
        self.len_T = 100  # [um], Tissue layer width
        self.len_VW = 10  # [um], Vessel wall width (from paper)
        self.R = self.int_r + self.len_EC + self.len_T + self.len_VW

        self.wss = 1.5  # [Pa], Wall shear stress
        self.wss_ref = 2.4  # [Pa], Reference wall shear stress

class Oxygen:
    def __init__(self):
        self.d_coeff = 2800  # [um^2/s], Diffusion coefficient for O2
        self.Q_max_vw = 5  # [uM/s], Max O2 consumption rate at vascular wall
        self.Q_max_t = 50  # [uM/s], Max O2 consumption rate at tissue
        self.Km_eNOS = 4.7  # [Torr]
        self.P = 70  # [Torr], P_O2 in blood lumen

class NitricOxide:
    def __init__(self, param):
        self.d_coeff = 3300  # [um^2/s], Diffusion coefficient for NO
        self.q_ref = 50  # [uM/s], Reference/control NO production rate
        self.C_ref = 27e-3  # [uM], Reference NO concentration
        self.R_max = param.wss / param.wss_ref * self.q_ref
