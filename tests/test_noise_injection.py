"""Summary."""
import numpy as np
import pytest

import smiter
from smiter.noise_functions import (
    GaussNoiseInjector,
    UniformNoiseInjector,
    JamssNoiseInjector,
)
from smiter.synthetic_mzml import Scan

# np.random.seed(1312)

# make it a test class
# integration test, test GaussNoiseInjector methods by themselves!
def test_GaussNoiseInjector_intensity_noise_ms1():
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": i_array,
            "ms_level": 1,
        }
    )
    noise_injector = GaussNoiseInjector(variance=0.2)
    scan = noise_injector.inject_noise(scan)
    # breakpoint()
    # assert abs((scan.i - 3e6) < 3e6 * 0.05).all()
    assert (scan.i - i_array < i_array * 0.2).all()


# make it a test class
# integration test, test GaussNoiseInjector methods by themselves!
def test_GaussNoiseInjector_intensity_noise_ms2():
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": i_array,
            "ms_level": 2,
        }
    )
    noise_injector = GaussNoiseInjector()
    scan = noise_injector.inject_noise(scan, dropout=0.0)
    assert (scan.i - i_array < i_array * 0.2).all()


def test_GaussNoiseInjector_intensity_noise_ms2_drop_all():
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": i_array,
            "ms_level": 2,
        }
    )
    noise_injector = GaussNoiseInjector()
    scan = noise_injector.inject_noise(scan, dropout=1.0)
    assert len(scan.mz) == 0
    assert len(scan.i) == 0


def test_GaussNoiseInjector_generate_intensity_noise():
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": i_array,
            "ms_level": 1,
        }
    )
    noise_injector = GaussNoiseInjector()
    noise = noise_injector._generate_intensity_noise(scan)
    assert (abs(noise) < i_array * 0.2).all()


def test_GaussNoiseInjector_generate_mz_noise():
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": i_array,
            "ms_level": 1,
        }
    )
    noise_injector = GaussNoiseInjector()
    noise = noise_injector._generate_mz_noise(scan)
    # breakpoint()
    assert (abs(noise) < scan.mz * 1e-5).all()


def test_gauss_mz_noise_ms1():
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": i_array,
            "ms_level": 1,
        }
    )
    noise_injector = GaussNoiseInjector()
    scan = noise_injector.inject_noise(scan)
    # breakpoint()
    assert (
        abs(scan.mz - np.array([100, 200, 300], dtype="float32"))
        < np.array([100, 200, 300], dtype="float32") * 1e-5
    ).all()


def test_uniform_intensity_noise_ms1():
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": i_array,
            "ms_level": 1,
        }
    )
    noise_injector = UniformNoiseInjector(ppm_noise=5e-6, intensity_noise=0.4)
    scan = noise_injector.inject_noise(scan)
    assert (abs(scan.i - i_array) < i_array * 0.4).all()


# make it a test class
# integration test, test UniformNoiseInjector methods by themselves!
def test_uniform_intensity_noise_ms2():
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": i_array,
            "ms_level": 2,
        }
    )
    noise_injector = UniformNoiseInjector(ppm_noise=5e-6, intensity_noise=0.4)
    scan = noise_injector.inject_noise(scan, dropout=0.0)
    # breakpoint()
    print("before", i_array)
    print("after", scan.i)
    print(abs(i_array - scan.i))
    print()
    assert (abs(scan.i - i_array) < i_array * 0.4).all()


def test_uniform_intensity_noise_ms2_drop_all():
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": i_array,
            "ms_level": 2,
        }
    )
    # breakpoint()
    noise_injector = UniformNoiseInjector(dropout=1.0)
    scan = noise_injector.inject_noise(scan)
    # breakpoint()
    assert len(scan.mz) == 0
    assert len(scan.i) == 0


def test_UniformNoiseInjector_generate_intensity_noise():
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": i_array,
            "ms_level": 1,
        }
    )
    noise_injector = UniformNoiseInjector(ppm_noise=5e-6, intensity_noise=0.4)
    noise = noise_injector._generate_intensity_noise(scan, intensity_noise=0.4)
    assert (abs(noise) < i_array * 0.4).all()


def test_UniformNoiseInjector_generate_mz_noise():
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": i_array,
            "ms_level": 1,
        }
    )
    noise_injector = UniformNoiseInjector(ppm_noise=5e-6, intensity_noise=0.4)
    noise = noise_injector._generate_mz_noise(scan, ppm_noise=5e-6)
    # breakpoint()
    assert (abs(noise) < scan.mz * 5e-6).all()


def test_UniformNoiseInjector_mz_noise_ms1():
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": i_array,
            "ms_level": 1,
        }
    )
    noise_injector = UniformNoiseInjector(ppm_noise=5e-6, intensity_noise=0.4)
    scan = noise_injector.inject_noise(scan)
    # breakpoint()
    assert (
        abs(scan.mz - np.array([100, 200, 300], dtype="float32"))
        < np.array([100, 200, 300], dtype="float32") * 5e-6
    ).all()


# def test_JamssNoiseInjector_mz_noise_ms1():
#    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
#    scan = Scan(
#        {
#            "mz": np.array([100, 200, 300], dtype="float32"),
#            "i": i_array,
#            "ms_level": 1,
#        }
#    )
#    noise_injector = JamssNoiseInjector(ppm_noise=5e-6, intensity_noise=0.4)
#    scan = noise_injector.inject_noise(scan)


def test_JamssNoiseInjector_intensity_noise_ms1():
    i_array = np.array([1e6, 2e6, 3e6], dtype="float32")
    scan = Scan(
        {
            "mz": np.array([100, 200, 300], dtype="float32"),
            "i": i_array,
            "ms_level": 1,
        }
    )
    noise_injector = JamssNoiseInjector(ppm_noise=5e-6, intensity_noise=0.4)
    scan = noise_injector.inject_noise(scan)
