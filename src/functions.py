import numpy as np
import scipy.constants as const

def fdd( E, mu, T): #fermi Dirac function
    if T == 0:
        f = np.heaviside(-(E-mu), 1)
    else:
        f = 1/(1+np.exp((E-mu)/(const.k*T/const.e)))
    return f


def dynesdos(E, Gamma, Delta): #dynes function
    dos = np.real((np.abs(E+1j*Gamma))/np.sqrt((E+1j*Gamma)**2-Delta**2))
    return np.abs(dos)

def dynesConvolute(V,E_int,conductance,delta,T,gamma):
    curr = []
    for Vp in V:
        currp = np.trapz((conductance)*dynesdos(E_int-Vp,gamma,delta)*(fdd(E_int, Vp, T)-fdd(E_int,0, T)),x=E_int)
        curr.append(currp)
    return np.gradient(np.array(curr))