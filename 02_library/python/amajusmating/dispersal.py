"""
Functions related to calculating probabilities of dispersal under the generalised
normal distribution.

Tom Ellis
"""

import numpy as np
from scipy.special import gamma
import pandas as pd

def distance_matrix(x, y):
    """
    Euclidean distance between two arrays.

    Parameters
    ----------
    x, y N-dimensional arrays.

    Returns
    -------
    x by y dimensional array.
    """
    distance_matrix = (x[np.newaxis] - y[:, np.newaxis])**2
    distance_matrix = distance_matrix.sum(axis=2)
    distance_matrix = np.sqrt(distance_matrix).T
    return distance_matrix


def generalised_normal_PDF(x, scale, shape, gamma_b=None):
    """
    Calculate the PDF of the generalised normal distribution.

    Parameters
    ----------
    x: vector
        Vector of deviates from the mean.
    scale: float
        Scale parameter.
    shape: float
        Shape parameter
    gamma_b: float, optional
        To speed up calculations, values for Euler's gamma for 1/shape
        can be calculated ahead of time and included as a vector.
    """
    xv = np.copy(x)
    if gamma_b:
        return (shape/(2 * scale * gamma_b ))      * np.exp(-(xv/scale)**shape)
    else:
        return (shape/(2 * scale * gamma(1.0/shape) )) * np.exp(-(xv/scale)**shape)

def dispersal_GND(x, scale, shape, w):
    """
    Calculate a probability that each candidate is a sire
    assuming assuming he is either drawn at random form the
    population, or from a generalised normal function of his
    distance from each mother. The relative contribution of the
    two distributions is controlled by mixture parameter w.

    Parameters
    ----------
    x: 2-d array
        Matrix of deviates from the mean.
    scale: float
        Scale parameter.
    shape: float
        Shape parameter
    w: float between 0 and 1.
        The proportion of probability mass assigned to the
        generalised normal function.
    """
    prob_GND = generalised_normal_PDF(x, scale, shape)
    prob_GND = prob_GND / prob_GND.sum(axis=1)[:, np.newaxis]

    prob_drawn = (prob_GND * w) + ((1-w) / x.shape[1])
    prob_drawn = prob_drawn / prob_drawn.sum(1, keepdims=True)
    prob_drawn = np.log(prob_drawn)

    return prob_drawn

def stdev_GND(scale, shape):
    """
    Calculate the 2nd moment of the genealised normal distribution.

    Parameters
    ----------
    scale: float
        Scale parameter.
    shape: float
        Shape parameter

    Returns
    -------
    Float.
    """
    return (scale * gamma(2.0/shape)) / (1.0/shape)

def grid_interpolation_GND(x, shape_vals, scale_vals):
    """
    Get ML estimates for the shape and scale parameters of the GND
    by grid interpolation.
    
    Parameters
    ----------
    x: 1-d array
        Vector of distances.
    shape_vals: 1-d array
        Vector of shape parameters.
    scale_vals: 1-d array
        Vector of scale parameters.
        
    Returns
    -------
    Dataframe showing ML vales for shape and scale, and the associated
    log likelihood.  
    
    """
    #likmat = np.empty([len(shape_vals), len(scale_vals)])
    likmat = []
    for shape in shape_vals:
        for scale in scale_vals:
            pr_dist = generalised_normal_PDF(x, scale = scale, shape = shape)
            pr_dist = np.log(pr_dist).sum()
            #likmat[shape, scale] = np.log(pr_dist).sum()
            likmat = likmat + [[shape, scale, pr_dist]]

    likmat = pd.DataFrame(likmat, columns=["shape", "scale", "loglik"])
    
    return likmat.loc[likmat['loglik'] == likmat['loglik'].max()]
