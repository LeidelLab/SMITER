from tempfile import NamedTemporaryFile

import pytest

import pymzml
from scipy.stats import kstest, normaltest
from smiter.synthetic_mzml import write_mzml


def simple_fragmentation_function(mol):
    """Test frag function.

    Args:
        mol (str): molecule to fragment

    Returns:
        list: list of fragment mzs
    """
    return [200]


def test_write_mzml():
    """Write a mzML without spectra readable by pymzML."""
    file = NamedTemporaryFile("wb")
    molecules = []
    peak_props = {}

    mzml_path = write_mzml(file, molecules, simple_fragmentation_function, peak_props)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 0


def test_write_inosine_flat_mzml():
    file = NamedTemporaryFile("wb")
    molecules = ["+C(10)H(12)N(4)O(5)"]
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": None,
        }
    }
    mzml_path = write_mzml(file, molecules, simple_fragmentation_function, peak_props)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 1001


def test_write_inosine_gauss_mzml():
    file = NamedTemporaryFile("wb")
    molecules = ["+C(10)H(12)N(4)O(5)"]
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 0.1},  # 10% of peak width
        }
    }
    mzml_path = write_mzml(file, molecules, simple_fragmentation_function, peak_props)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 1001
    intensities = []
    for spec in reader:
        if spec.ms_level == 1:
            intensities.append(spec.i[0])
    # assert peaks are gauss distributed
    _, p = normaltest(intensities)
    assert p < 5e-4


def test_write_inosine_gamma_mzml():
    file = NamedTemporaryFile("wb")
    molecules = ["+C(10)H(12)N(4)O(5)"]
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gamma",
            "peak_params": {"a": 3, "scale": 20},  # 10% of peak width,
        }
    }
    mzml_path = write_mzml(file, molecules, simple_fragmentation_function, peak_props)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 1001
    intensities = []
    for spec in reader:
        if spec.ms_level == 1:
            intensities.append(spec.i[0])
    # assert peaks are gauss distributed
    t = kstest(intensities, "gamma", args=(3, 0, 20))
    assert t.pvalue < 5e-4
