import numpy as np
import PathLog #FIXME this needs to come from somewhere


# FIXME: There is a definite weirdness in the definition of the run functions
def run_case_1(method, E, dE, x0, k,interp, tol, maxnit, precon_scheme, path_traverse,
        fixed_ends, verbose=method, precon, precon_prep = precon_scheme, direction = path_traverse):

    # initialize variables
    #FIXME: the julia code says x = x0.x, but I have no idea what the equivalent of . is in python

    nit = 0
    numdE = 0
    numE = 0

    log = PathLog()
    if verbose >= 2:
        pass #FIXME put something here

    file = []
    if verbose >= 4:
        pass #FIXME put something here

    xout, log, alpha = odesolve(solver(method))
