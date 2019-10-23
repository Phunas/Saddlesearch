
import numpy as np

def rayleigh(v, x, len, E0, E, P):
    w = v/(np.sqrt(np.dot(v,np.matmul(P,v))))
    return 2.0 * (E(x+len/2*w) - 2.0 * E0 + E(x-len/2*v)) / (len^2)

lambda localmerit(x, x0, v0, len, g0, λ0, E) = (0.5 * ( E(x+len/2*v0) + E(x-len/2*v0) ) - 2.0 * dot(v0, g0) * dot(v0, x-x0)- λ0 * np.dot(v0, x-x0)^2  )


