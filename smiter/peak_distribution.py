#!/usr/bin/env python3
"""Distribution funtion for chromo peaks.

Attributes:
    distributions (dict): mapping distribution name to distribution function
"""
import math
from typing import Callable, Dict

from loguru import logger
from scipy.stats import gamma


def gauss_dist(x: float, sigma: float = 1, mu: float = 0):
    """Calc Gauss distribution.

    Args:
        x (float): x
        sigma (float, optional): standard deviation
        mu (float, optional): mean

    Returns:
        float: y
    """
    # logger.debug(f"Sigma: {sigma}")
    # print(f'Sigma: {sigma}')
    return (
        (
            1
            / (sigma * math.sqrt(2 * math.pi))
            * pow(math.e, (-0.5 * pow(((x - mu) / sigma), 2)))
        )
        / 0.3989422804014327
        * sigma
    )


def gauss_tail(
    x: float,
    mu: float,
    sigma: float,
    scan_start_time: float,
    h: float = 1,
    t: float = 0.2,
    f: float = 0.01,
) -> float:
    # h = 1
    ## should be a parameter
    # t = 0.2
    # f = 0.01
    sigma = t * (x - scan_start_time) + f  # / (x + 0.00001)
    i = h * math.e ** (-0.5 * ((x - mu) / sigma) ** 2)
    return i


def gamma_dist(x: float, a: float = 5, scale: float = 0.33):
    """Calc gamma distribution.

    Args:
        x (float): Description
        a (float, optional): Description
        scale (float, optional): Description

    Returns:
        float: y
    """
    return gamma.pdf(x, a=a, scale=scale)


distributions = {
    "gauss": gauss_dist,
    "gamma": gamma_dist,
    "gauss_tail": gauss_tail,
}  # type: Dict[str, Callable]
