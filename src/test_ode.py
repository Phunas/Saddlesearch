import numpy as np

from ode_old_new_in_progress import odesolve_r12


def x_2(x, nit):
    x_np1 = x
    A = np.array([[1,0],[0,0.01]])
    x_2 = np.dot(A, x_np1)
    return -x_2

print(odesolve_r12(x_2, np.array([5, 5]), np.array([]), 1))
