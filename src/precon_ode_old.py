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
"""


def odesolve_r12(f, X0, atoms, h=None, verbose=1, fmax=1e-6, maxtol=1e3,
             steps=100, rtol=1e-1, 
                 C1 = 1e-2, C2=2.0, hmin=1e-10, extrapolate=3, callback=None, precon=None):

    X = X0
    X_out = [] #Create an array to store the values of X
    
    Fn = f(X) #Get the Forces
    if precon is not None:

        print(Fn)

        precon.make_precon(atoms)
        Fn = precon.solve(Fn)
        print(Fn)
 
    Rn = np.linalg.norm(Fn, np.inf) #Get the residual: This could be a different norm. Note always goes to zero because forces
    
    X_out.append(X) #appending current values of X
    log = []
    log.append([0, X, Rn]) #appending current residual
    
    if Rn <= fmax:
        print(f"ODE12r terminates succesfully after 0 iterations")#Forces are already small
        return X_out, log, h

    if Rn >= maxtol:
        print(f"SADLESEARCH: Residual {Rn} is too large at itteration 0") #Forces are too big
        return X_out, log, h

    # computation of the initial step
    r = np.linalg.norm(Fn, np.inf) #pick the biggest force
    if h is None:
        h = 0.5 * rtol**0.5 / r #Chose a stepsize based on that force
        h = max(h, hmin) #Make sure the step size is not too big

    for nit in range(1, steps):
        #Redistribute
        Xnew = X + h * Fn #Pick a new position
        Fnew = f(Xnew) # Calculate the new forces at this position

       
        if precon is not None:
            atoms_new = atoms.copy()
            atoms_new.set_positions(Xnew.reshape(len(atoms_new), 3))
            precon.make_precon(atoms_new)
            Fnew = precon.solve(Fnew)


        Rnew = np.linalg.norm(Fnew, np.inf) #Find the new residual forces

        e = 0.5 * h * (Fnew - Fn) #Estimate the area under the foces curve
       
        err = np.linalg.norm(e, np.inf) # Come up with an error based on this area

        #This deceides whether or not to acccept the new residual 
        if Rnew <= Rn*(1-C1*h) or Rnew <= (Rn*C2 and err <= rtol ):
            accept = True

        else:
            accept = False
            conditions = (Rnew <= Rn * (1 - C1 * h), Rnew <= Rn * C2, err <= rtol ) # THIS ALSO SEEMS POTENTIALLY WRONG


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
        if np.isnan(h_ls) or h_ls < hmin: #This rejects the increment if it is too small or the extrapolation scheme misbehaves
            h_ls = np.inf

        # This pickes a separate increment
        h_err = h * 0.5 * np.sqrt(rtol/err)

        #We incremnet the system
        if accept:
            X = Xnew
            Fn = Fnew

            Rn = Rnew
            if callback is not None:
                callback(X)
            
            X_out = np.append(X_out, X) #Store X
            log = np.append(log, Rn) 

            #We check the residuals again
            if Rn <= fmax:
                if verbose >= 1:
                    print(f"SADDLESEARCH: terminates succesfully after {nit} iterations.")
                X_out = np.append(X_out, X)
                log = np.append(log, Rn) 
                return X_out, log, h
            if Rn >= maxtol:
                print(f"SADLESEARCH: Residual {Rn} is too large at itteration number {nit}")
        
                X_out = np.append(X_out, X) #Store X
                log = np.append(log, Rn) 
                return X_out, log, h

            # Compute a new step size. This is based on the extrapolation and some other heuristics
            h = max(0.25 * h, min(4*h, h_err, h_ls))
            # Log step-size analytic results

        else:
        # Compute a new step size.
            h = max(0.1 * h, min(0.25*h, h_err, h_ls)) #This also computes a new step size if the old one is not allowed 
        # error message if step size is too small
        if abs(h) <= hmin:
            print(f'ODE12r Step size {h} too small at nit = {nit}')
            return X_out, log, h


    #Logging:
    if verbose >= 1:
        print(f'ODE12r terminates unuccesfully after {steps} iterations.')

    return X_out, log, h

from ase.optimize.sciopt import SciPyOptimizer

from ase.optimize.precon import Exp, C1, Pfrommer

class ODE12rOptimizer_precon(SciPyOptimizer):
    def call_fmin(self, fmax, steps):

        X_out, log, h = odesolve_r12(lambda x: -self.fprime(x),
                                     self.x0(),
                                     self.atoms,
                                     fmax=fmax,
                                     steps=steps,
                                     callback=self.callback,
                                     precon = Exp())
        
