import numpy as np

from ode_old_new_in_progress import odesolve_r12


def x_2(x, P, targ):
    x = float(x[-1])
    x_2 = x*x
    Rn = x_2 - targ
    return x_2, Rn, -1, np.dot 



print(odesolve_r12(x_2, np.array([5]), np.array([]), 1))
