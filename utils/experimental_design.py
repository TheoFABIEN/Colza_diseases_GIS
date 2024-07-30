import numpy as np 
from scipy.stats import qmc 


def scaled_latin_hypercube(param_ranges, sample_size):
    """
    Generate a scaled Latin Hypercube sample within specified parameter ranges.

    Args:
        param_ranges (numpy.array): A 2D array where each row specifies the 
        minimum and maximum values for a parameter.
        sample_size (int): The number of samples to generate.

    Returns:
        numpy.array: A 2D array of shape (sample_size, len(param_ranges)) 
        containing the scaled Latin Hypercube samples.

    Example:
        param_ranges = np.array([[0, 1], [10, 20], [100, 200], [-50, 50]])
        sample_size = 10
        samples = scaled_latin_hypercube(param_ranges, sample_size)
    """
    sampler = qmc.LatinHypercube(d = len(param_ranges))
    samples = sampler.random(n = sample_size)

    scaled_samples = qmc.scale(samples, param_ranges[:,0], param_ranges[:,1])

    return scaled_samples

