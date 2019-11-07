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


def odesolve_r12(f, X0, log, P, g = lambda P, X: X, precon_prep= lambda P, X: X, file=None, verbose=1, tol=1e-4, maxtol=1e3,
             maxint=100, method='ODE', rtol=1e-1, threshold=1,
             C1 = 1e-2, C2=2.0, h=None, hmin=1e-10,MaxF=1e3, Extrapolate=3):
    if verbose >= 2:
        pass #FIXME put something here

    if verbose >= 4 and file != None:
        pass # FIXME put someting here

    X = X0
    P = precon_prep(P, X)
    #print(X)
    X_out = []

    numdE = 0
    numE = 0

    X = g(X,P)
    
    P = precon_prep(P, X)
    Fn, Rn, ndE, _ = f(X,P,0)

    numdE += ndE

    np.append(X0, X)

    np.append(log, numE)
    np.append(log, numdE)
    np.append(log, Rn)

    if verbose >= 2:
        pass #FIXME put something here
    if verbose >= 4 and file != None:
        pass # FIXME put someting here

    if Rn <= tol:
        if verbose >= 1:
            print("SADDLESEARCH: {} terminates succesfully after {} iterations.".format(method, nit))
        if verbose >=4:
            Pass #FIXME Something Goes here

        return X_out, log, h


        return X_out, log, h

    if Rn >= maxtol:
        print("SADLESEARCH: Residual {} is too large at itteration number {}".format(Rn, nit))
        if verbose >= 4:
            pass #FIXME something Goes here.
        np.append(X_out, X) #Store X
        np.append(log, numE)
        np.append(log, numdE)
        np.append(log, Rn)
        return X_out, log, h

#computation of the initail step
#    print(Fn)
#    print(X[-1])
#    print(np.absolute(X[-1]))
#    print(np.max(np.absolute(X[-1])))
#    print(np.divide(Fn, np.max(np.absolute(X[-1]))))
    #r = np.linalg.norm(np.array([np.divide(Fn, np.max(np.absolute(X[-1]), threshold))]), inf) 
    r = np.abs(np.divide(Fn, np.max(np.absolute(X[-1]))))
    if h == None:
        h = 0.5 * rtol**0.5 /r
        h = max(h, hmin)

    for nit in range(1, maxint):
        #Redistribute
        Xnew = g(X+h*Fn, P)  # that wayNone it implicitly becomes part
                                # of `f`
                                # but it seems to make the evolution slower; need more testing!
        #return force
        Pnew = precon_prep(P, Xnew)
        Fnew, Rnew, ndE, dot_P = f(Xnew, Pnew, nit)

        numdE +=ndE

        # error estimation

        e = 0.5 * h * (Fnew - Fn)
       
        print('blah')
        #print(e)
        #np.array([np.abs(X), np.abs(Xnew)])
        #blah = max(np.array([np.abs(X), np.abs(Xnew)]))
        #print(e, np.max(np.array([np.abs(X), np.abs(Xnew)]))) 
        #print(np.divide(e,max(np.max(np.array([np.abs(float(X)), np.abs(float(Xnew))])))))
        #err = np.linalg.norm(np.divide(e,np.max(np.max([np.abs(X), np.abs(Xnew)]), threshold)), 'inf') # THIS LINE FEELS WRONG
        #err = np.linalg.norm(np.divide(e, np.max(np.ndarray([np.abs(float(X)), np.abs(float(Xnew))]))), 'inf')
        #blah = np.array([np.abs(X), np.abs(Xnew)])
        blah = float(max(np.array([np.abs(X), np.abs(Xnew)])))
        blah = np.divide(e, float(max(np.array([np.abs(X), np.abs(Xnew)]))))
        blah = np.divide(np.array([e]), max(np.array([np.abs(X), np.abs(Xnew)])))
        print('blah', blah)
        err = np.linalg.norm(np.array(blah), np.inf)
        print(err)
        #1/0
        # accept step if residual is sufficient decreased

        if Rnew <= Rn*(1-C1*h) or Rnew <= (Rn*C2 and err <= rtol ):
            accept = True
        else:
            accept = False
            conditions = (Rnew <= Rn * (1 - C1 * h), Rnew <= Rn * C2, err <= rtol ) # THIS ALSO SEEMS POTENTIALLY WRONG

        # whether we accept or reject this step, we now need a good guess for
        # the next step-size, from a line-search-like construction
        y = Fn - Fnew
        if Extrapolate == 1:       # F(xn + h Fn) ⋅ Fn ~ 0
            h_ls = h*dot_P(Fn, y)/(dot_P(y,y))
        elif Extrapolate == 2:   # F(Xn + h Fn) ⋅ F{n+1} ~ 0
            h_ls = h * dot_P(Fn, Fnew) / (dot_P(Fn, y) + 1e-10)
        elif Extrapolate == 3:   # min | F(Xn + h Fn) |
            h_ls = h * dot_P(Fn, y) / (dot_P(y, y) + 1e-10)
        else:
            print('SADDLESEARCH: invalid extrapolate parameter')
            if verbose >= 4 and file != None:
                pass #FIXME include VERBOSE options
            raise NameError('invalid extrapolate parameter')
        if np.isnan(h_ls) or h_ls < hmin:
            h_ls = np.inf

        h_err = h * 0.5 * np.sqrt(rtol/err)


        if accept == True:
            X = Xnew
            Fn = Fnew
            Rn = Rnew
            P = Pnew


        np.append(X_out, X) #Store X
        np.append(log, numE)
        np.append(log, numdE)
        np.append(log, Rn) 
        if verbose >=2:
            pass #FIXME include VERBOSE options
        if verbose >= 4:
            pass #FIXME include VERBOSE options

        if Rn <= tol:
            if verbose >= 1:
                print("SADDLESEARCH: {} terminates succesfully after {} iterations.".format(method, nit))
            if verbose >= 4 and file!= None:
                pass # FIXME put someting here
            return X_out, log, h

        if Rn >= maxtol:
            print("SADLESEARCH: Residual {} is too large at itteration number {}".format(Rn, nit))
            if verbose >= 4 and file!= None:
                pass # FIXME put someting here
        
            np.append(X_out, X) #Store X
            np.append(log, numE)
            np.append(log, numdE)
            np.append(log, Rn) 
            return X_out, log, h

        # Compute a new step size.
        h = max(0.25 * h, min(4*h, h_err, h_ls))
        # Log step-size analytic resukts
        if verbose >= 3:
            pass # FIXME put someting here
        if verbose >= 4 and file!= None:
            pass # FIXME put someting here

    else:
        # Compute a new step size.
        h = max(0.1 * h, min(4*h, h_err, h_ls))
        if verbose >= 3:
            pass # FIXME put someting here
        if verbose >= 4 and file != None:
            pass # FIXME put someting here

    # error message if step size is too small
    if abs(h) <= hmin:
        print('SADLESEARCH: Step size {} too small at nit = {}'.format(h, nit))
        if verbose >= 4 and file != None:
            pass # FIXME put someting here
        return X_out, log, h


    #Logging:
    if verbose>=1:
        print('SADDLESEARCH {} terminates unuccesfully after {} iterations.'.format(method, maxnit))

    if verbose>=4 and file != None:
        pass #FIXME: Insert something here

    return X_out, log, h
####### STOP WHEN YOU HAVE REACHED 383
