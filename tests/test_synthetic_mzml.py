from tempfile import NamedTemporaryFile

import numpy as np
import pymzml
from scipy.signal import find_peaks
from scipy.stats import kstest, normaltest

import pytest
import smiter
from smiter.fragmentation_functions import (
    AbstractFragmentor,
    NucleosideFragmentor,
    PeptideFragmentor,
)
from smiter.noise_functions import GaussNoiseInjector, UniformNoiseInjector
from smiter.synthetic_mzml import write_mzml


class TestFragmentor(AbstractFragmentor):
    def __init__(self):
        print("Fragmentor goes chop chop chop chop chop ...")

    def fragment(self, mol):
        return np.array([(200, 1e5)])


fragmentor = TestFragmentor()
noise_injector = GaussNoiseInjector(variance=0.05)


def test_write_mzml():
    """Write a mzML without spectra readable by pymzML."""
    file = NamedTemporaryFile("wb")
    peak_props = {}
    mzml_params = {
        "gradient_length": 0,
    }
    mzml_path = write_mzml(file, peak_props, fragmentor, noise_injector, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 0


def test_write_empty_mzml():
    """Write a mzML without spectra readable by pymzML."""
    file = NamedTemporaryFile("wb")
    peak_props = {}
    mzml_params = {
        "gradient_length": 5,
    }
    mzml_path = write_mzml(file, peak_props, fragmentor, noise_injector, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 166


def test_write_inosine_flat_mzml():
    file = NamedTemporaryFile("wb")
    peak_props = {
        "inosine": {
            "chemical_formula": "+C(10)H(12)N(4)O(5)",
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
    mzml_path = write_mzml(file, peak_props, fragmentor, noise_injector, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 999


def test_write_inosine_gauss_mzml():
    file = NamedTemporaryFile("wb")
    peak_props = {
        "inosine": {
            "chemical_formula": "+C(10)H(12)N(4)O(5)",
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
    mzml_path = write_mzml(file, peak_props, fragmentor, noise_injector, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 999
    intensities = []
    for spec in reader:
        if spec.ms_level == 1:
            if len(spec.i) > 0:
                intensities.append(spec.i[0])
    _, p = normaltest(intensities)
    print(intensities)
    # assert p < 5e-4


def test_write_mzml_get_TIC():
    file = NamedTemporaryFile("wb")
    peak_props = {
        "inosine": {
            "chemical_formula": "+C(10)H(12)N(4)O(5)",
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
    noise_injector = GaussNoiseInjector(variance=0.0)
    mzml_path = write_mzml(file, peak_props, fragmentor, noise_injector, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 999
    reader = pymzml.run.Reader(mzml_path)
    tics = []
    for spec in reader:
        if spec.ms_level == 1:
            tics.append(sum(spec.i))
    tic = reader["TIC"]
    assert tic.peaks()[:, 1] == pytest.approx(np.array(tics, dtype="float32"))


def test_write_inosine_gamma_mzml():
    file = NamedTemporaryFile("wb")
    peak_props = {
        "inosine": {
            "chemical_formula": "+C(10)H(12)N(4)O(5)",
            "trivial_name": "inosine",
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gamma",
            "peak_scaling_factor": 1e6,
            "peak_params": {"a": 3, "scale": 20},  # 10% of peak width,
        }
    }
    mzml_params = {
        "gradient_length": 30,
    }
    mzml_path = write_mzml(file, peak_props, fragmentor, noise_injector, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 999
    intensities = []
    for spec in reader:
        if spec.ms_level == 1:
            if len(spec.i) > 0:
                intensities.append(spec.i[0])
    ## what is a reasonable p-value cutoff here?
    # assert t.pvalue < 1e-100


def test_write_inosine_adenosine_gauss_mzml():
    file = NamedTemporaryFile("wb")
    peak_props = {
        "inosine": {
            "chemical_formula": "+C(10)H(12)N(4)O(5)",
            "charge": 2,
            "trivial_name": "inosine",
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
        },
        "+C(10)H(13)N(5)O(4)": {
            "chemical_formula": "+C(10)H(13)N(5)O(4)",
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
    mzml_path = write_mzml(file, peak_props, fragmentor, noise_injector, mzml_params)
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
    # assert p < 5e-4
    _, p = normaltest(adeno_intensities)
    # assert p < 5e-4


def test_write_inosine_adenosine_gauss_shift_mzml():
    file = NamedTemporaryFile("wb")
    peak_props = {
        "inosine": {
            "chemical_formula": "+C(10)H(12)N(4)O(5)",
            "charge": 2,
            "trivial_name": "inosine",
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
            "peak_scaling_factor": 1e6,
        },
        "+C(10)H(13)N(5)O(4)": {
            "chemical_formula": "+C(10)H(13)N(5)O(4)",
            "trivial_name": "adenosine",
            "charge": 2,
            "scan_start_time": 15,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
            "peak_scaling_factor": 1e6,
        },
    }
    mzml_params = {
        "gradient_length": 45,
    }
    noise_injector = GaussNoiseInjector(variance=0.0)
    # noise_injector = UniformformNoiseInjector()
    mzml_path = write_mzml(file, peak_props, fragmentor, noise_injector, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 1499

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
    # assert rt max diff is about 15
    m_i = np.argmax(ino_intensities)
    m_a = np.argmax(adeno_intensities)
    mean_i_rt = ino_rt[m_i]
    mean_a_rt = adeno_rt[m_a]
    assert 14 < (mean_a_rt - mean_i_rt) < 16


def test_write_inosine_proper_fragments_mzml():
    file = NamedTemporaryFile("wb")
    peak_props = {
        "inosine": {
            "chemical_formula": "+C(10)H(12)N(4)O(5)",
            "charge": 2,
            "trivial_name": "inosine",
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
            "peak_scaling_factor": 1e6,
        },
    }
    mzml_params = {
        "gradient_length": 30,
    }

    nucl_fragmentor = NucleosideFragmentor()
    mzml_path = write_mzml(
        file, peak_props, nucl_fragmentor, noise_injector, mzml_params
    )
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
            # 9 out of 10 specs are empty, only collect the one
            if len(spec.mz) > 0:
                ino_fragments.append(spec.mz)
            pass

    _, p = normaltest(ino_intensities)
    # assert p < 5e-4
    expected_frags = np.array([137.0457872316])
    # Check ino fragments are correct
    from pprint import pprint

    # pprint(ino_fragments)
    for frag_list in ino_fragments:
        print(frag_list)
        assert len(frag_list) == len(expected_frags)
        sorted_frags = np.sort(frag_list)
        assert abs(sorted_frags - expected_frags) < 0.001


@pytest.mark.slow()
def test_write_peptide_gauss_mzml():
    file = NamedTemporaryFile("wb")
    peak_props = {
        "ELVISLIVES": {
            "chemical_formula": "ELVISLIVES",
            "trivial_name": "ELVISLIVES",
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width
        },
        "ELVISLIVSE": {
            "charge": 2,
            "trivial_name": "ELVISLIVSE",
            "chemical_formula": "ELVISLIVSE",
            "scan_start_time": 15,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
        },
    }
    mzml_params = {
        "gradient_length": 45,
    }
    fragmentor = PeptideFragmentor()
    mzml_path = write_mzml(file, peak_props, fragmentor, noise_injector, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 1499
    intensities = []
    for spec in reader:
        if spec.ms_level == 1:
            if len(spec.i) > 0:
                intensities.append(sum(spec.i))
    _, p = normaltest(intensities)
    # assert p < 5e-4


def test_write_2_mols_same_cc():
    file = NamedTemporaryFile("wb")
    peak_props = {
        "uridine": {
            "charge": 2,
            "chemical_formula": "+C(9)H(11)N(2)O(6)",
            "trivial_name": "uridine",
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 1},  # 10% of peak width,
            "peak_scaling_factor": 0.5 * 1e6,
        },
        "pseudouridine": {
            "chemical_formula": "+C(9)H(11)N(2)O(6)",
            "trivial_name": "pseudouridine",
            "charge": 2,
            "scan_start_time": 15,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 1},  # 10% of peak width,
            "peak_scaling_factor": 1e6,
        },
    }
    mzml_params = {
        "gradient_length": 45,
        "min_intensity": 10,
    }
    mzml_params = {
        "gradient_length": 45,
    }
    fragmentor = NucleosideFragmentor()
    noise_injector = GaussNoiseInjector(variance=0.0)
    mzml_path = write_mzml(file, peak_props, fragmentor, noise_injector, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    # assert reader.get_spectrum_count() == 1499
    intensities = []
    for spec in reader:
        if spec.ms_level == 1:
            if len(spec.i) > 0:
                intensities.append(sum(spec.i))
    peaks, _ = find_peaks(intensities)
    assert len(peaks) == 2


def test_rescale_intensity():
    i = 100
    rt = 15
    molecule = ""
    peak_props = {
        "uridine": {
            "charge": 2,
            "chemical_formula": "+C(9)H(11)N(2)O(6)",
            "trivial_name": "uridine",
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 1},  # 10% of peak width,
            "peak_scaling_factor": 0.5,
        }
    }
    iso_lib = {}

    i_rescaled = smiter.synthetic_mzml.rescale_intensity(
        i, rt, "uridine", peak_props, iso_lib
    )
    # should be max_i @ half peak width, with half peak width = 15 and max_i = 50 (100 * 0.5)
    assert round(i_rescaled, 2) == 50


def test_generate_scans():
    pass


def test_generate_molecule_isotopologue_lib():
    pass


def test_write_scans():
    pass
