import numpy as np

"""
`ODE12r`: adaptive ODE solver, uses 1st and 2nd order approximations to estimate local error and find a new step length
### Parameters:
* `rtol` : relative tolerance
* `threshold` : threshold for error estimate
* `C1` : sufficient contraction parameter
* `C2` : residual growth control (Inf means there is no control)
* `h` : step size, if nothing is passed, an estimate is used based on ODE12
* `hmin` : minimal allowed step size
* `maxF` : terminate if |Fn| > maxF * |F0|
* `extrapolate` : extrapolation style (3 seems the most robust)
"""


def odesolve_r12(f, X0, log, P, h=0.2, g = lambda P, X: X, precon_prep= lambda P, X: X, file=None, verbose=1, tol=1e-6, maxtol=1e3,
             maxint=100, method='ODE', rtol=1e-1, threshold=1,
             C1 = 1e-2, C2=2.0, hmin=1e-10,MaxF=1e3, Extrapolate=3):

    X = X0
    X_out = []
    
    Fn = f(X,0)
    print("Fn", Fn)
    Rn = np.linalg.norm(Fn, np.inf)
    
    X_out = np.append(X_out, X)
    log = np.append(log, Rn)
    
    if Rn <= tol:
        print("Done")
        #print("SADDLESEARCH: {} terminates succesfully after {} iterations.".format(method, nit))
        return X_out, log, h

    if Rn >= maxtol:
        print("SADLESEARCH: Residual {} is too large at itteration number {}".format(Rn, nit))
        #np.append(X_out, X) #Store X
        #np.append(log, Rn)
        return X_out, log, h

    # computation of the initial step
    r = np.linalg.norm(Fn, np.inf)
    if h == None:
        h = 0.5 * rtol**0.5 /r
        h = max(h, hmin)

    for nit in range(1, maxint):
#        print(nit)
        #Redistribute
        Xnew = X + h * Fn
        print(X, Xnew)
        Fnew = f(Xnew, nit)
        Rnew = np.linalg.norm(Fnew, np.inf)

        e = 0.5 * h * (Fnew - Fn)
       
        err = np.linalg.norm(e, np.inf)

        if Rnew <= Rn*(1-C1*h) or Rnew <= (Rn*C2 and err <= rtol ):
            accept = True
            print("accept = True")
        else:
            print("accept = False")
            accept = False
            conditions = (Rnew <= Rn * (1 - C1 * h), Rnew <= Rn * C2, err <= rtol ) # THIS ALSO SEEMS POTENTIALLY WRONG

        y = Fn - Fnew
        if Extrapolate == 1:       # F(xn + h Fn) ⋅ Fn ~ 0
            h_ls = h*np.dot(Fn, y)/(np.dot(y,y))
        elif Extrapolate == 2:   # F(Xn + h Fn) ⋅ F{n+1} ~ 0
            h_ls = h * np.dot(Fn, Fnew) / (np.dot(Fn, y) + 1e-10)
        elif Extrapolate == 3:   # min | F(Xn + h Fn) |
            h_ls = h * np.dot(Fn, y) / (np.dot(y, y) + 1e-10)
#            print("Fn", Fn)
#            print("y", y)             
#            print("h_ls", h_ls)
#            1/0
        else:
            print('SADDLESEARCH: invalid extrapolate parameter')
            raise NameError('invalid extrapolate parameter')
#        if np.isnan(h_ls.any()) or np.all(h_ls < hmin):
        if np.isnan(h_ls) or h_ls < hmin:
            h_ls = np.inf

        h_err = h * 0.5 * np.sqrt(rtol/err)


        if accept == True:
            X = Xnew
            Fn = Fnew
            Rn = Rnew


            X_out = np.append(X_out, X) #Store X
            log = np.append(log, Rn) 

            if Rn <= tol:
                if verbose >= 1:
                    print("SADDLESEARCH: {} terminates succesfully after {} iterations.".format(method, nit))
                X_out = np.append(X_out, X)
                log = np.append(log, Rn) 
                return X_out, log, h
            if Rn >= maxtol:
                print("SADLESEARCH: Residual {} is too large at itteration number {}".format(Rn, nit))
        
                X_out = np.append(X_out, X) #Store X
                log = np.append(log, Rn) 
                return X_out, log, h

            # Compute a new step size.
            h = max(0.25 * h, min(4*h, h_err, h_ls))
            # Log step-size analytic resukts

        else:
            print("false")
        # Compute a new step size.
            h = max(0.1 * h, min(0.25*h, h_err, h_ls))
            print(h)
            print(h_err)
            print(h_ls)
        # error message if step size is too small
        if abs(h) <= hmin:
            print('SADLESEARCH: Step size {} too small at nit = {}'.format(h, nit))
            return X_out, log, h


    #Logging:
    if verbose>=1:
        print("Done2")
        print('SADDLESEARCH {} terminates unuccesfully after {} iterations.'.format(method, maxint))

    return X_out, log, h
####### STOP WHEN YOU HAVE REACHED 383
