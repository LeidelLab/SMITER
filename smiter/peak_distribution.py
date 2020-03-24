#!/usr/bin/env python3
"""Distribution funtion for chromo peaks.

Attributes:
    distributions (dict): mapping distribution name to distribution function
"""
import math

from scipy.stats import gamma


def gauss_dist(x, sigma=1, mu=0):
    """Calc Gauss distribution.

    Args:
        x (str): x
        sigma (float, optional): standard deviation
        mu (float, optional): mean

    Returns:
        float: y
    """
    return (
        1
        / (sigma * math.sqrt(2 * math.pi))
        * pow(math.e, (-0.5 * pow(((x - mu) / sigma), 2)))
    )


def gamma_dist(x, a=5, scale=0.33):
    """Calc gamma distribution.

    Args:
        x (float): Description
        a (float, optional): Description
        scale (float, optional): Description

    Returns:
        float: y
    """
    return gamma.pdf(x, a=a, scale=scale)


def flat(x):
    """Summary.

    Args:
        x (TYPE): Description

    Returns:
        TYPE: Description
    """
    return x


distributions = {
    "gauss": gauss_dist,
    "gamma": gamma_dist,
    "flat": flat,
}
