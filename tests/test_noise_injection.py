"""Summary."""
import numpy as np
import pytest

from smiter.noise_functions import TestNoiseInjector
from smiter.synthetic_mzml import Scan

np.random.seed(1312)

# make it a test class
# integration test, test TestNoiseInjector methods by themselves!
def test_simple_intensity_noise_ms1():
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": np.array([1e6, 2e6, 3e6], dtype="float32"),
            "ms_level": 1,
        }
    )
    noise_injector = TestNoiseInjector()
    scan = noise_injector.inject_noise(scan)
    assert abs((scan.i - 3e6) < 3e6 * 0.05).all()


# make it a test class
# integration test, test TestNoiseInjector methods by themselves!
def test_simple_intensity_noise_ms2():
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": np.array([1e6, 2e6, 3e6], dtype="float32"),
            "ms_level": 2,
        }
    )
    noise_injector = TestNoiseInjector()
    scan = noise_injector.inject_noise(scan)
    assert abs((scan.i - 3e6) < 3e6 * 0.05).all()


def test_TestNoiseInjector_generate_intensity_noise():
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": np.array([1e6, 2e6, 3e6], dtype="float32"),
            "ms_level": 1,
        }
    )
    noise_injector = TestNoiseInjector()
    noise = noise_injector._generate_intensity_noise(scan)
    assert (abs(noise) < np.array([1e6, 2e6, 3e6], dtype="float32") * 0.05).all()


def test_TestNoiseInjector_generate_mz_noise():
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": np.array([1e6, 2e6, 3e6], dtype="float32"),
            "ms_level": 1,
        }
    )
    noise_injector = TestNoiseInjector()
    noise = noise_injector._generate_mz_noise(scan)
    # breakpoint()
    assert (abs(noise) < scan.mz * 1e-5).all()


def test_simple_mz_noise_ms1():
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": np.array([1e6, 2e6, 3e6], dtype="float32"),
            "ms_level": 1,
        }
    )
    noise_injector = TestNoiseInjector()
    scan = noise_injector.inject_noise(scan)
    # breakpoint()
    assert (
        abs(scan.mz - np.array([100, 200, 300], dtype="float32"))
        < np.array([100, 200, 300], dtype="float32") * 1e-5
    ).all()
