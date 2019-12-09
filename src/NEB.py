import numpy as np


def run_case_1(method, x0, k,interp, tol, maxnit, precon_scheme, path_traverse,
        fixed_ends, verbose=method, precon, precon_prep = precon_scheme, direction = path_traverse):

    #HOW DO WE INITIALIZE VARIABLES JULIA LINE 7

    xout, alpha = odesolve(solver(method), lambda X: forces(typex0,X,direction(NI, nit), k, interp, fixed_ends)), tol = tol, maxtol = maxtol, maxnit = maxnit, method = "$(typeof(method))")



### STARTING BACK AT LINE 80 FOR THE CONTINUATION, WILL COME BACK TO IT

def FORCES(path_type, X, dE, 
           direcion, k, interp, fixed_ends)


    if interp == 1:
        # central finite differences
        dxds = [[np.zeros(x[1])], [0.5*(x[i+1]-x[i-1]) for i=2:N-1], [np.zeros(x[1])]]

    elif interp > 1:
        # splines
        dxds = [dist(precon_scheme, P, x, i) for i=1:length(x)-1]
        param = [0; [sum(ds[1:i]) for i in 1:length(ds)]]
        param /= param[end]
        param[end] = 1
        d2xds2 = parametrise(dxxs, x, ds, parametrisation = param)
        k *= (1/(N**2))
    else:
        print('SADDLESEARCH: invalid `interpolate` parameter')

    dxds ./= point_norm(precon_scheme, P, dxds)
    dxds[1] = zeros(dxds[1])
    dxds[-1] =  = zeros(dxds[1])

    # elastic interactions between adjacent images
    Fk = elastic_force(precon_scheme, P, k*N*M, dxds, d2xds2)

    # potential gradient
    dE0_temp = []

    if fixed_ends != True:
        dE0_temp = [dE(x[i]) for i in direction]
        cost = N
    else:
        dE0_temp = [[zeros(x[1])]; [dE(x[i]) for i in direction[2:end-1]]; [zeros(x[1])]]
        cost = N-2

    dE0 = [dE0_temp[i] for i in direction]


    # projecting out tangent term of potential gradient
    dE0_perp = proj_grad(precon_scheme, Pm de0, dxds)

    # collecting force term
    F = forcing(precon_scheme, precon, dE0_perp-Fk)

    # residual error
    res = maxres(precon_scheme, P, dE0_perp)

    return F, res, cost, (SOMETHING GOES HERE), P, convert(path_type, Y))

def jacobian(precon, path_type, X, dE, ddE, k)

    x = convert(path_type, X)

    #preconditioner
    Np = size(precon, 1)
    P_1 = lambda i: precon[np.mod(i-1, Np) +1, 1]
    P_2 = lambda i,j: precon{np.mod(i-1,Np) +1, mod(j-1,Np)+1]

    #FIXME TWO LINES I CANNOT PROPERLY TRANSLATE hessian = ddE.precon; hessian_prep! = ddE.precon_prep!
                                                 #hessian = hessian_prep!(hessian, x)


    H_1(i) = hessian[np.mod(i-1,Np)+1, 1]
    H_2(i,j) = hessian[np.mod(i-1,Np)+1, np.mod(j-1,Np)+1]

    N = length(x)
    M = length(x[1])
    O = zeros(M, M)
    J = np.ndarray.fill(O,(N, N))

    [J[n,n-1] = dF_nm(x,n,dE,P) + dS

