import numpy as np

"""
`ODE12r`: adaptive ODE solver, uses 1st and 2nd order approximations to estimate local error and find a new step length
### Parameters:
* `rtol` : relative tolerance
* `C1` : sufficient contraction parameter
* `C2` : residual growth control (Inf means there is no control)
* `h` : step size, if nothing is passed, an estimate is used based on ODE12
* `hmin` : minimal allowed step size
* `maxF` : terminate if |Fn| > maxF * |F0|
* `extrapolate` : extrapolation style (3 seems the most robust)
* `res_norm` : The norm that we pick for calculating the residual
"""

def odesolve_r12(f, X0, h=None, verbose=1, fmax=1e-6, maxtol=1e3,
                 steps=100, rtol=1e-1, res_norm=np.inf, error_norm=np.inf, step_norm=np.inf,
                 C1 = 1e-2, C2=2.0, hmin=1e-10, extrapolate=3, callback=None):
    X = X0
    X_out = [] #Create an array to store the values of X
    
    Fn = f(X) #Get the Forces
    Rn = np.linalg.norm(Fn, res_norm) #Get the residual.
    
    X_out.append([X]) 
    quick_log = []
    log = []
    quick_log.append([0, Rn])
    log.append([0, X, Rn]) 
    
    if Rn <= fmax:
        print("ODE12r terminates succesfully at initial step.") #System already converged
        return X_out, log, h

    if Rn >= maxtol:
        print(f"ODE12r: Residual {Rn} is too large at inintial step") #Forces are too big
        return X_out, quicklog, h

    # computation of the initial step
    r = np.linalg.norm(Fn, step_norm) #pick the biggest force
    if h is None:
        h = 0.5 * rtol**0.5 / r #Chose a stepsize based on that force
        h = max(h, hmin) #Make sure the step size is not too big

    for nit in range(1, steps):
        #Redistribute
        Xnew = X + h * Fn
        Fnew = f(Xnew)
        Rnew = np.linalg.norm(Fnew, res_norm)

        e = 0.5 * h * (Fnew - Fn) #Estimate the area under the foces curve
       
        err = np.linalg.norm(e, error_norm) # Come up with an error based on this area

        #This deceides whether or not to acccept the new residual 
        if Rnew <= Rn*(1-C1*h) or Rnew <= (Rn*C2 and err <= rtol ):
            accept = True
        else:
            accept = False
            conditions = (Rnew <= Rn * (1 - C1 * h), Rnew <= Rn * C2, err <= rtol)

        #This decides on an extrapolation scheme for the system, to pick a new increment.
        y = Fn - Fnew
        if extrapolate == 1:       # F(xn + h Fn)
            h_ls = h*np.dot(Fn, y)/(np.dot(y,y))
        elif extrapolate == 2:   # F(Xn + h Fn)
            h_ls = h * np.dot(Fn, Fnew) / (np.dot(Fn, y) + 1e-10)
        elif extrapolate == 3:   # min | F(Xn + h Fn) |
            h_ls = h * np.dot(Fn, y) / (np.dot(y, y) + 1e-10)
        else:
            raise ValueError(f'invalid extrapolate parameter: {extrapolate}')
        if np.isnan(h_ls) or h_ls < hmin
            h_ls = np.inf

        h_err = h * 0.5 * np.sqrt(rtol/err)

        if accept:
            X = Xnew
            Fn = Fnew
            Rn = Rnew
            if callback is not None:
                callback(X)
            
            X_out.append([X]) #Store X
            log.append([nit, Rn, X]) 
            quicklog.append([nit, Rn])

            #We check the residuals again
            if Rn <= fmax:
                if verbose >= 1:
                    print(f"ODE12r: terminates succesfully after {nit} iterations.")
                X_out.append([X])
                log.append([nit, Rn, X])
                quicklog.append([nit, Rn])
                return X_out, quicklog, h
            if Rn >= maxtol:
                print(f"ODE12r: Residual {Rn} is too large at itteration number {nit}")
        
                X_out.append([X])
                log.append([nit, Rn, X])
                quicklog.append([nit, Rn])
                return X_out, quicklog, h

            h = max(0.25 * h, min(4*h, h_err, h_ls))

        else:
            h = max(0.1 * h, min(0.25*h, h_err, h_ls)) 
        if abs(h) <= hmin:
            print(f'SADLESEARCH: Step size {h} too small at nit = {nit}')
            return X_out, quicklog, h

    #Logging:
    if verbose >= 1:
        print(f'ODE12r terminates unuccesfully after {nit} iterations.')

    return X_out, quicklog, h

from ase.optimize.sciopt import SciPyOptimizer

class ODE12rOptimizer(SciPyOptimizer):
    def call_fmin(self, fmax, steps):
        X_out, log, h = odesolve_r12(lambda x: -self.fprime(x),
                                     self.x0(),
                                     fmax=fmax,
                                     steps=steps,
                                     callback=self.callback)
        
