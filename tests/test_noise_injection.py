"""Summary."""
import numpy as np
import pytest

from smiter.noise_functions import TestNoiseInjector
from smiter.synthetic_mzml import Scan

# make it a test class
# integration test, test TestNoiseInjector methods by themselves!
def test_simple_intensity_noise_ms1():
    scan = Scan({
        "mz": np.array([100, 200, 300]),
        "i":  np.array([1e6, 2e6, 3e6]),
        "ms_level": 1,
        })
    noise_injector = TestNoiseInjector()
    scan = noise_injector.inject_noise(scan)
    abs(scan.i - 3e6) < 3e6 * 0.05


# make it a test class
# integration test, test TestNoiseInjector methods by themselves!
def test_simple_intensity_noise_ms2():
    scan = Scan({
        "mz": np.array([100, 200, 300]),
        "i":  np.array([1e6, 2e6, 3e6]),
        "ms_level": 2,
        })
    noise_injector = TestNoiseInjector()
    scan = noise_injector.inject_noise(scan)
    abs(scan.i - 3e6) < 3e6 * 0.05


def test_TestNoiseInjector_generate_intensity_noise():
    scan = Scan({
        "mz": np.array([100, 200, 300]),
        "i":  np.array([1e6, 2e6, 3e6]),
        "ms_level": 1,
    })
    noise_injector = TestNoiseInjector()
    noise = noise_injector._generate_intensity_noise(scan)
    # breakpoint()
    assert noise == pytest.approx(np.array([ -22509.97615172, -171353.03222904,   98737.02411275], dtype='float32'))


def test_TestNoiseInjector_generate_mz_noise():
    scan = Scan({
        "mz": np.array([100, 200, 300]),
        "i":  np.array([1e6, 2e6, 3e6]),
        "ms_level": 1,
    })
    noise_injector = TestNoiseInjector()
    noise = noise_injector._generate_mz_noise(scan)
    assert noise is None
