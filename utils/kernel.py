from scipy.special import gamma
import numpy as np 

def dispersion_kernel(x, y, mu = 250, beta2 = 0.2, normalize_kernel = True):
    """
    Computes a 2D exponential power dispersal kernel based on the given x and 
    y coordinates. If `normalize_kernel` is True, the kernel is normalized by 
    the sum of its values along the y-axis.

    Args:
        x (numpy.array): 1D array of x-coordinates.
        y (numpy.array): 1D array of y-coordinates.
        mu (float): Mean dispersion distance of the kernel. Default is 250.
        beta2 (float): Shape parameter for the kernel. Default is 0.2.
        normalize_kernel (bool): If True (default), the kernel is normalized. 

    Returns:
        numpy.array: The 2D dispersion kernel.
    """
    x = x[:, None]
    beta1 = mu*gamma(2/beta2)/gamma(3/beta2)
    res = beta2/(2*np.pi*beta1*gamma(2/beta2))
    res*= np.exp(-(np.sqrt(x*x + y*y)/beta1)**beta2)

    if normalize_kernel:
        w = np.apply_along_axis(func1d = np.trapz, axis = 1, arr = res)
        w[np.isnan(w)] = 0
        res/= sum(w)
        return res

    return res
    
    
    
    
    
    
    
    
    
