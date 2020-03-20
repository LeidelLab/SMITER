#!/usr/bin/env python3
import math

from scipy.stats import gamma


def gauss_dist(x, sigma=1, mu=0):
    return (
        1
        / (sigma * math.sqrt(2 * math.pi))
        * pow(math.e, (-0.5 * pow(((x - mu) / sigma), 2)))
    )


def gamma_dist(x, a=5, scale=0.33):
    return gamma.pdf(x, a=a, scale=scale)


def flat(x):
    return x


distributions = {
    "gauss": gauss_dist,
    "gamma": gamma_dist,
    "flat": flat,
}
