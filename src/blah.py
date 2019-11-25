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


@with_kw type ODE12r
   rtol::Float64 = 1e-1
   threshold::Float64 = 1.0
   C1::Float64 = 1e-2
   C2::Float64 = 2.0
   h = nothing
   hmin::Float64 = 1e-10
   maxF::Float64 = 1e3
   extrapolate::Int = 3
end

function odesolve(solver::ODE12r, f, X0::Vector{Float64}, log::IterationLog;
                  file = nothing,
                  verbose = 1,
                  g=(X, P)->X, tol=1e-4, maxtol=1e3, maxnit=100,
                  P = I, precon_prep! = (P, X) -> P,
                  method = "ODE" )

   X = copy(X0)
   Xout = [];

   # initialise variables
   Fn, Rn, ndE, _ = f(X, P, 0)

   push!(Xout, X)
   push!(log, numE, numdE, Rn)

   # logging
   if Rn <= tol
      if verbose >= 1
         println("SADDLESEARCH: $method terminates succesfully after $(nit) iterations.")
      end
      return Xout, log, h
   end
   if Rn >= maxtol
      warn("SADDLESEARCH: Residual $Rn is too large at nit = $nit.");
      push!(Xout, X) # store X
      push!(log, typemax(Int64), typemax(Int64), Rn) # residual, store history
      return Xout, log, h
   end

   # computation of the initial step
   r = norm(Fn ./ max.(abs.(X), threshold), Inf) + realmin(Float64)
   if h == nothing
      h = 0.5 * rtol^(1/2) / r
      h = max(h, hmin)
   end

   for nit = 1:maxnit
      Xnew = g(X + h * Fn, P)
      Fnew, Rnew, ndE, dot_P = f(Xnew, Pnew, nit)
      # error estimation
      e = 0.5 * h * (Fnew - Fn)
      err = norm(e ./ max.(maximum([abs.(X) abs.(Xnew)],2), threshold), Inf) + realmin(Float64)

      # accept step if residual is sufficient decreased
      if (   ( Rnew <= Rn * (1 - C1 * h) )         # contraction
          || ( Rnew <= Rn * C2 && err <= rtol ) )  # moderate growth + error control
         accept = true
      else
         accept = false
         conditions = (Rnew <= Rn * (1 - C1 * h), Rnew <= Rn * C2, err <= rtol )
      end

      # whether we accept or reject this step, we now need a good guess for
      # the next step-size, from a line-search-like construction
      y = Fn - Fnew
      if extrapolate == 1       # F(xn + h Fn) ⋅ Fn ~ 0
         h_ls = h * dot_P(Fn, Fn) / dot_P(Fn, y)
      elseif extrapolate == 2   # F(Xn + h Fn) ⋅ F{n+1} ~ 0
         h_ls = h * dot_P(Fn, Fnew) / (dot_P(Fn, y) + 1e-10)
      elseif extrapolate == 3   # min | F(Xn + h Fn) |
         h_ls = h * dot_P(Fn, y) / (dot_P(y, y) + 1e-10)
      else
         @printf("SADDLESEARCH: invalid `extrapolate` parameter")
         error("SADDLESEARCH: invalid `extrapolate` parameter")
      end
      if isnan(h_ls) || (h_ls < hmin)
         h_ls = Inf
      end
      # or from the error estimate
      h_err = h * 0.5 * sqrt(rtol/err)

      if accept
         X, Fn, Rn, P  = Xnew, Fnew, Rnew, Pnew

         push!(Xout, X) # store X
         push!(log, numE, numdE, Rn) # residual, store history

         if Rn <= tol
            if verbose >= 1
               println("SADDLESEARCH: $(method) terminates succesfully after $(nit) iterations.")
            end
            return Xout, log, h
         end

         if Rn >= maxtol
            warn("SADDLESEARCH: Residual $Rn is too large at nit = $nit.");
            push!(Xout, X) # store X
            push!(log, typemax(Int64), typemax(Int64), Rn) # residual, store history
            return Xout, log, h
         end

         # Compute a new step size.
         h = max(0.25 * h, min(4*h, h_err, h_ls))
         # log step-size analytic results
      else
         # compute new step size
         h = max(0.1 * h, min(0.25 * h, h_err, h_ls))
         # log step-size analytic results
      end

      # error message if step size is too small
      if abs(h) <= hmin
         warn("SADDLESEARCH: Step size $h too small at nit = $nit.");
         return Xout, log, h

   # logging
   if verbose >= 1
      println("SADDLESEARCH: $(method) terminated unsuccesfully after $(maxnit) iterations.")

   return Xout, log, h
end
