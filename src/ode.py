import numpy as np



class ODE(Calculator):

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

    def __init__()


