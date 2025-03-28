import autograd.numpy as np  # Use autograd's numpy instead of regular numpy
from autograd import grad

# Constants
Cp = 1005.7  # from Davies-Jones function, was 1005.
Cpv = 1870.0  # J/kg/K
Cw = 4190.0
L0 = 2.501e6  # J/kg

C = 273.15  # K
Rd = 287.04
Rv = 461.5
RdoRv = Rd / Rv


def LvK(TempK):
    """Latent heat of water vapor
    
    Args:
        TempK: Temperature in Kelvin
        
    Returns:
        Latent heat in J/kg
    """
    return L0 + (Cpv - Cw) * (TempK - 273.0)


def es(T, p=1e5, TK=None, P=None):
    """Saturation vapor pressure based on Wexler's formula,
    with enhancement factor for moist air rather than water vapor.
    The enhancement factor requires a pressure.
    
    Args:
        T: Temperature in degrees C (ignored if TK is provided)
        p: Pressure in Pa (default 1e5 Pa)
        TK: Temperature in Kelvin (optional)
        P: Pressure in hPa (optional)
        
    Returns:
        Saturation vapor pressure in Pa
        
    From A. L. Buck 1981: JAM, 20, 1527-1532.
    """
    if P is not None:
        p_hPa = P
    else:
        p_hPa = p * 1e-2
        
    if TK is not None:
        T = TK - C
        
    # Use autograd.numpy.exp instead of numpy.exp
    esat = 1e2 * 6.1121 * (1.0007 + 3.46e-8 * p_hPa) * np.exp((17.502 * T) / (240.97 + T))
    return esat  # in Pa


def qs(p, T):
    """Saturation specific humidity based on Wexler's formula for es
    with enhancement factor.
    
    Args:
        p: Pressure in Pa
        T: Temperature in degrees C
        
    Returns:
        Saturation specific humidity in kg/kg
        
    From A. L. Buck 1981: JAM, 20, 1527-1532.
    """
    esat = es(T, p)
    return RdoRv * esat / (p + (RdoRv - 1) * esat)


def dqsdT(p, T):
    """Derivative of qs with respect to T at p,T by autodiff
    
    Args:
        p: Pressure in Pa
        T: Temperature in degrees C
        
    Returns:
        Derivative of saturation specific humidity with respect to temperature
    """
    # Use autograd to compute the derivative
    return grad(lambda t: qs(p, t))(T)


def updatex(f, x, fhat=0):
    """General single Newton iteration to update x toward f(x) = fhat for a univariate function f
    
    Args:
        f: Function to find root of
        x: Current x value
        fhat: Target function value (default 0)
        
    Returns:
        Updated x value
    """
    f_prime = grad(f)
    return x + (fhat - f(x)) / f_prime(x)


def Twet_autodiff(T, q, p, niter=2):
    """Wet bulb temperature using Newton's method for target specific humidity q.
    Uses automatic differentiation.
    
    Args:
        T: Temperature in Kelvin
        q: Specific humidity in kg/kg
        p: Pressure in Pa
        niter: Number of Newton iterations (default 2)
        
    Returns:
        Wet bulb temperature in Kelvin
    """
    def f(t):
        return (t - T) + LvK((T + t) / 2) / Cp * (qs(p, t - C) - q)
    
    t = T
    for i in range(niter):
        t = updatex(f, t, 0)
    
    return t