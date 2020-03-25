from tempfile import NamedTemporaryFile

import numpy as np
import pymzml
import pytest
from scipy.stats import kstest, normaltest

from smiter.fragmentation_functions import AbstractFragmentor, NucleosideFragmentor
from smiter.synthetic_mzml import write_mzml


def simple_fragmentation_function(mol):
    """Test frag function.

    Args:
        mol (str): molecule to fragment

    Returns:
        list: list of fragment mzs
    """
    return [200]


class TestFragmentor(AbstractFragmentor):
    def __init__(self):
        print("Fragmentor goes chop chop chop chop chop ...")

    def fragment(self, mol):
        return np.array([(200, 1e5)])


fragmentor = TestFragmentor()


def test_write_mzml():
    """Write a mzML without spectra readable by pymzML."""
    file = NamedTemporaryFile("wb")
    molecules = []
    peak_props = {}
    mzml_params = {
        "gradient_length": 0,
    }
    mzml_path = write_mzml(file, molecules, fragmentor, peak_props, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 0


def test_write_inosine_flat_mzml():
    file = NamedTemporaryFile("wb")
    molecules = ["+C(10)H(12)N(4)O(5)"]
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "trivial_name": "inosine",
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": None,
            "peak_params": {},
        }
    }
    mzml_params = {
        "gradient_length": 30,
    }
    mzml_path = write_mzml(file, molecules, fragmentor, peak_props, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 999


def test_write_inosine_gauss_mzml():
    file = NamedTemporaryFile("wb")
    molecules = ["+C(10)H(12)N(4)O(5)"]
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "trivial_name": "inosine",
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width
        }
    }
    mzml_params = {
        "gradient_length": 30,
    }
    mzml_path = write_mzml(file, molecules, fragmentor, peak_props, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 999
    intensities = []
    for spec in reader:
        if spec.ms_level == 1:
            intensities.append(spec.i[0])
    _, p = normaltest(intensities)
    print(intensities)
    assert p < 5e-4


def test_write_inosine_gamma_mzml():
    file = NamedTemporaryFile("wb")
    molecules = ["+C(10)H(12)N(4)O(5)"]
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "trivial_name": "inosine",
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gamma",
            "peak_params": {"a": 3, "scale": 20},  # 10% of peak width,
        }
    }
    mzml_params = {
        "gradient_length": 30,
    }
    mzml_path = write_mzml(file, molecules, fragmentor, peak_props, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 999
    intensities = []
    for spec in reader:
        if spec.ms_level == 1:
            intensities.append(spec.i[0])
    t = kstest(intensities, "gamma", args=(3, 0, 20))
    print(intensities)
    ## what is a reasonable p-value cutoff here?
    assert t.pvalue < 1e-100


def test_write_inosine_adenosine_gauss_mzml():
    file = NamedTemporaryFile("wb")
    molecules = ["+C(10)H(12)N(4)O(5)", "+C(10)O(4)N(5)H(13)"]
    # molecules = ["+C(10)O(4)N(5)H(13)"]
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "charge": 2,
            "trivial_name": "inosine",
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
        },
        "+C(10)H(13)N(5)O(4)": {
            "trivial_name": "adenosine",
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
        },
    }
    mzml_params = {
        "gradient_length": 30,
    }
    mzml_path = write_mzml(file, molecules, fragmentor, peak_props, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 999

    ino_rt = []
    ino_intensities = []
    adeno_rt = []
    adeno_intensities = []

    ino_mono = 269.0880
    adeno_mono = 268.1040

    for spec in reader:
        if spec.ms_level == 1:
            ino_i = spec.i[abs(spec.mz - ino_mono) < 0.001]
            adeno_i = spec.i[abs(spec.mz - adeno_mono) < 0.001]
            if len(adeno_i) > 0:
                adeno_intensities.append(adeno_i[0])
                adeno_rt.append(spec.scan_time[0])
            if len(ino_i) > 0:
                ino_intensities.append(ino_i[0])
                ino_rt.append(spec.scan_time[0])
    _, p = normaltest(ino_intensities)
    assert p < 5e-4
    _, p = normaltest(adeno_intensities)
    assert p < 5e-4


def test_write_inosine_adenosine_gauss_shift_mzml():
    file = NamedTemporaryFile("wb")
    molecules = ["+C(10)H(12)N(4)O(5)", "+C(10)O(4)N(5)H(13)"]
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "charge": 2,
            "trivial_name": "inosine",
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
        },
        "+C(10)H(13)N(5)O(4)": {
            "trivial_name": "adenosine",
            "charge": 2,
            "scan_start_time": 15,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
        },
    }
    mzml_params = {
        "gradient_length": 45,
    }
    mzml_path = write_mzml(file, molecules, fragmentor, peak_props, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 1507

    ino_rt = []
    ino_intensities = []
    adeno_rt = []
    adeno_intensities = []

    ino_mono = 269.0880
    adeno_mono = 268.1040

    for spec in reader:
        if spec.ms_level == 1:
            ino_i = spec.i[abs(spec.mz - ino_mono) < 0.001]
            adeno_i = spec.i[abs(spec.mz - adeno_mono) < 0.001]
            if len(adeno_i) > 0:
                adeno_intensities.append(adeno_i[0])
                adeno_rt.append(spec.scan_time[0])
            if len(ino_i) > 0:
                ino_intensities.append(ino_i[0])
                ino_rt.append(spec.scan_time[0])
    _, p = normaltest(ino_intensities)
    assert p < 5e-4
    _, p = normaltest(adeno_intensities)
    assert p < 5e-4

    # assert rt max diff is about 15
    m_i = np.argmax(ino_intensities)
    m_a = np.argmax(adeno_intensities)
    mean_i_rt = ino_rt[m_i]
    mean_a_rt = adeno_rt[m_a]
    assert 14 < (mean_a_rt - mean_i_rt) < 16


def test_write_inosine_proper_fragments_mzml():
    file = NamedTemporaryFile("wb")
    molecules = ["+C(10)H(12)N(4)O(5)"]
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "charge": 2,
            "trivial_name": "inosine",
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
        },
    }
    mzml_params = {
        "gradient_length": 30,
    }

    nucl_fragmentor = NucleosideFragmentor()

    mzml_path = write_mzml(file, molecules, nucl_fragmentor, peak_props, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 999

    ino_rt = []
    ino_intensities = []
    ino_fragments = []

    ino_mono = 269.0880

    for spec in reader:
        if spec.ms_level == 1:
            ino_i = spec.i[abs(spec.mz - ino_mono) < 0.001]
            if len(ino_i) > 0:
                ino_intensities.append(ino_i[0])
                ino_rt.append(spec.scan_time[0])
        if spec.ms_level == 2:
            # check if inosine precursor
            # if true, add fragments to ino_fragments
            ino_fragments.append(spec.mz)
            pass

    _, p = normaltest(ino_intensities)
    assert p < 5e-4

    expected_frags = np.array([137.0457872316])
    # Check ino fragments are correct
    from pprint import pprint
    # pprint(ino_fragments)
    for frag_list in ino_fragments:
        assert len(frag_list) == len(expected_frags)
        # breakpoint()
        sorted_frags = np.sort(frag_list)
        assert abs(sorted_frags - expected_frags) < 0.001
