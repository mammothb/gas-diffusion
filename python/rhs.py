import numpy as np

def make_o2_rhs(param, gas_o2, gas_no, u, v, r_coeff_o2, nr_12, r_23, r_34,
                r_45):
    # Extra reaction terms for the RHS
    R_NO_23 = gas_no.R_max * v[r_23] / (v[r_23] + gas_o2.Km_eNOS)
    app_Km_34 = param.Km * (1 + u[r_34] / gas_no.C_ref)
    R_O2_34 = gas_o2.Q_max_vw * v[r_34] / (v[r_34] + app_Km_34)
    app_Km_45 = param.Km * (1 + u[r_45] / gas_no.C_ref)
    R_O2_45 = gas_o2.Q_max_t * v[r_45] / (v[r_45] + app_Km_45)
    # RHS by layer
    C_12 = np.zeros(nr_12)
    C_23 = r_coeff_o2 * R_NO_23
    C_34 = r_coeff_o2 * R_O2_34
    C_45 = r_coeff_o2 * R_O2_45

    # Overall RHS
    return np.concatenate(([gas_o2.P], C_12, C_23, C_34, C_45))

def make_no_rhs(param, gas_o2, gas_no, u, v, r_coeff_no, r_01, nr_12, r_23,
                r_34, r_45, lambda_core):

    # Extra reaction terms for the RHS
    R_NO_01 = lambda_core * u[r_01]
    R_NO_23 = -gas_no.R_max * v[r_23] / (v[r_23] + gas_o2.Km_eNOS)
    R_NO_34 = param.lambda_vw * u[r_34]
    R_NO_45 = param.lambda_t * u[r_45]
    # RHS by layer
    B_01 = r_coeff_no * R_NO_01
    B_12 = np.zeros(nr_12)
    B_23 = r_coeff_no * R_NO_23
    B_34 = r_coeff_no * R_NO_34
    B_45 = r_coeff_no * R_NO_45

    # Overall RHS
    B = np.concatenate((B_01, B_12, B_23, B_34, B_45))
    B[0] /= 2
    return B
