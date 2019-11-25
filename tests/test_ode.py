import numpy as np

from ode_old_new_in_progress import odesolve_r12

def steepest_descent(f, x, h, steps=100, fmax=1e-6):
    log = []
    for i in range(steps):
        F = f(x, i)
        R = np.linalg.norm(F, np.inf)
        print(f'{i:8d} {x[0]:12.6f} {R:14.6e}')
        if R < fmax:
            return x, log
        x[:] = x + h*F
        log.append((x, R))
    raise RuntimeError('did not converge')
        

def x_2(x, nit):
    x_2 = 2*x
    return -x_2

x0 = np.array([5.0])
x = x0.copy()
x_sd, log_sd = steepest_descent(x_2, x, h=0.114)

x = x0.copy()
x_ode, log_ode, h = odesolve_r12(x_2, x, np.array([]), 1)
