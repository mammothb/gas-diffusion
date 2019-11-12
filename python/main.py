import math
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

import parameters
from rhs import make_no_rhs, make_o2_rhs

np.set_printoptions(threshold=np.inf)

def max_absolute_percentage_error(x, y):
    return np.max(np.abs((x - y) / np.maximum(np.finfo(float).eps, y)))

def main():
    param = parameters.Parameters()
    gas_no = parameters.NitricOxide(param)
    gas_o2 = parameters.Oxygen()

    # Simulation parameters
    h = 0.5  # space step
    omega = 1.96  # factor for seccessive overrelaxation method
    tolerance = 1e-6  # tolerance for relative error for Gauss-Siedel
    max_cfl = 5  # [um], Maximum CFL width to test for

    # Various compartments and domain space
    nr = round(param.R / h) + 1
    r = np.linspace(0, param.R, nr)

    # Big matrix for all solutions
    u_ans = np.zeros((nr, max_cfl))
    v_ans = np.zeros((nr, max_cfl))

    # Set up coefficient for R terms
    r_coeff_no = h**2 / gas_no.d_coeff
    r_coeff_o2 = h**2 / gas_o2.d_coeff / param.alpha

    # a-terms for LHS matrix
    a = h / 2 / r[1 : -1]
    # Create LHS for NO
    NO_shape = (nr, nr)
    D_NO = diags(-2 * np.ones(nr), 0, NO_shape)
    L_NO = diags(np.append(1 - a.T, 2), -1, NO_shape)
    U_NO = diags(np.append(2, 1 + a.T), 1, NO_shape)
    M_NO = L_NO + D_NO / omega
    N_NO = D_NO / omega - D_NO - U_NO

    for i in range(max_cfl):
        cfl = i + 1
        lambda_core = param.lambda_b / 2 * (1 + (param.int_r /
                                                 (param.int_r + cfl))**2)
        is_unsteady = True

        # Various compartments and domain space
        ind_r1 = int((param.int_r - cfl) / h)    # end of RBC core
        ind_r2 = int(param.int_r / h)            # end of vessel interior
        ind_r3 = int(ind_r2 + param.len_EC / h)  # end of EC layer
        ind_r4 = int(ind_r3 + param.len_VW / h)  # end of VW layer

        # Number of nodes in selected compartments
        nr_01 = ind_r1       # number of nodes for r0 <= r < r1
        nr_12 = ind_r2 - ind_r1  # number of nodes for r1 <= r < r2
        nr_15 = nr - nr_01       # number of nodes for r1 <= r <= r5

        # Index vectors for selected compartments for clarity
        r_01 = np.arange(0, ind_r1)       # r0 <= r < r1
        r_23 = np.arange(ind_r2, ind_r3)  # r2 <= r < r3
        r_34 = np.arange(ind_r3, ind_r4)  # r3 <= r < r4
        r_45 = np.arange(ind_r4, nr)      # r4 <= r <= r5

        # Initialization
        u = np.zeros(nr)

        v_01 = gas_o2.P * np.ones(nr_01)
        v = np.append(v_01, np.zeros(nr_15))

        # Display interface position
        print(f"CFL width: {cfl}")
        print("   r0    r1    r2    r3    r4    r5")
        print(f"{r[0]:5.1f} {r[ind_r1]:5.1f} {r[ind_r2]:5.1f} "
              f"{r[ind_r3]:5.1f} {r[ind_r4]:5.1f} {r[-1]:5.1f}")

        # Create LHS for O2
        O2_shape = (nr_15 + 1, nr_15 + 1)
        D_O2 = diags(np.append(1, -2 * np.ones(nr_15)), 0, O2_shape)
        L_O2 = diags(np.append(1 - a[ind_r1 - 1 :].T, 2), -1, O2_shape)
        U_O2 = diags(np.append(0, 1 + a[ind_r1 - 1 :].T), 1, O2_shape)
        M_O2 = L_O2 + D_O2 / omega
        N_O2 = D_O2 / omega - D_O2 - U_O2

        # Solve for an initial O2 profile
        G = make_o2_rhs(param, gas_o2, gas_no, u, v, r_coeff_o2, nr_12, r_23,
                        r_34, r_45)

        v_15 = spsolve(M_O2, N_O2 @ v[ind_r1 - 1 :] + G)
        v_new = np.concatenate((v_01, v_15[1 :]))

        while is_unsteady:
            # Solving for NO
            F = make_no_rhs(param, gas_o2, gas_no, u, v, r_coeff_no, r_01,
                            nr_12, r_23, r_34, r_45, lambda_core)
            u_new = spsolve(M_NO, N_NO @ u + F)
            # Solving for O2
            G = make_o2_rhs(param, gas_o2, gas_no, u, v, r_coeff_o2, nr_12,
                            r_23, r_34, r_45)
            v_15 = spsolve(M_O2, N_O2 @ v[ind_r1 - 1 :] + G)
            v_new = np.concatenate((v_01, v_15[1 :]))
            # Calculate relative errors
            u_err = max_absolute_percentage_error(u_new, u)
            v_err = max_absolute_percentage_error(v_new, v)
            # Output relative errors per iteration
            # print(u_err, v_err)
            if u_err < tolerance and v_err < tolerance:
                is_unsteady = False
            u = u_new
            v = v_new
        u_ans[:, i] = u
        v_ans[:, i] = v

    __, (axes_1, axes_2) = plt.subplots(2)
    for i in range(max_cfl):
        axes_1.plot(r, u_ans[:, i])
        axes_2.plot(r, v_ans[:, i])

    plt.show()

if __name__ == "__main__":
    start_time = time.time()
    main()
    print(time.time() - start_time)
