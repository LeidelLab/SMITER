from tempfile import NamedTemporaryFile

import numpy as np
import pymzml
import pytest
from scipy.signal import find_peaks
from scipy.stats import kstest, normaltest

import smiter
from smiter.fragmentation_functions import (
    AbstractFragmentor,
    NucleosideFragmentor,
    PeptideFragmentor,
)
from smiter.noise_functions import GaussNoiseInjector, UniformNoiseInjector
from smiter.synthetic_mzml import (
    generate_interval_tree,
    generate_molecule_isotopologue_lib,
    generate_scans,
    write_mzml,
)


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
    assert reader.get_spectrum_count() == 167


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
    assert reader.get_spectrum_count() == 1000


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
    assert reader.get_spectrum_count() == 1000
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
    assert reader.get_spectrum_count() == 1000
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
    assert reader.get_spectrum_count() == 1000
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
            "charge": 1,
            "trivial_name": "inosine",
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
            "peak_scaling_factor": 1e3,
        },
        "+C(10)H(13)N(5)O(4)": {
            "chemical_formula": "+C(10)H(13)N(5)O(4)",
            "trivial_name": "adenosine",
            "charge": 1,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
        },
    }
    mzml_params = {"gradient_length": 30, "min_intensity": 0}
    mzml_path = write_mzml(file, peak_props, fragmentor, noise_injector, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    assert reader.get_spectrum_count() == 1000

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
            "charge": 1,
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
            "charge": 1,
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
    assert reader.get_spectrum_count() == 1500

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
            "charge": 1,
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
    assert reader.get_spectrum_count() == 1000

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

    _, p = normaltest(ino_intensities)
    expected_frags = np.array([137.0457872316])

    for frag_list in ino_fragments:
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
    assert reader.get_spectrum_count() == 1500
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


def test_generate_interval_tree():
    peak_props = {
        "uridine": {
            "charge": 2,
            "chemical_formula": "+C(9)H(11)N(2)O(6)",
            "trivial_name": "uridine",
            "scan_start_time": 0,
            "peak_width": 5,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 1},  # 10% of peak width,
            "peak_scaling_factor": 1,
        },
        "pseudouridine": {
            "charge": 2,
            "chemical_formula": "+C(9)H(11)N(2)O(6)",
            "trivial_name": "uridine",
            "scan_start_time": 4,
            "peak_width": 5,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 1},  # 10% of peak width,
            "peak_scaling_factor": 1,
        },
    }
    tree = generate_interval_tree(peak_props)

    # only uridine
    tree_3 = tree[3]
    assert len(tree_3) == 1
    for entry in tree_3:
        assert entry.data == "uridine"
        assert entry.length() == 5

    # only pseudouridine
    tree_8 = tree[8]
    assert len(tree_3) == 1
    for entry in tree_8:
        assert entry.data == "pseudouridine"
        assert entry.length() == 5

    # uridine and pseudouridine
    tree_4 = tree[4]
    assert len(tree_4) == 2
    items = ["pseudouridine", "uridine"]
    for entry in tree_4:
        items.remove(entry.data)
        assert entry.length() == 5
    assert len(items) == 0


def test_generate_scans_simple():
    # generate 1 MS1 and 1 MS2 scan
    # peak width is as long as gradient and two times ms_rt_diff
    peak_props = {
        "uridine": {
            "charge": 2,
            "chemical_formula": "+C(9)H(11)N(2)O(6)",
            "trivial_name": "uridine",
            "scan_start_time": 0,
            "peak_width": 0.06,  # peak as long as ms_rt_diff and gradient length
            "peak_function": "gauss",
            "peak_params": {"sigma": 1},  # 10% of peak width,
            "peak_scaling_factor": 1,
        }
    }
    trivial_names = {val["chemical_formula"]: key for key, val in peak_props.items()}
    scans, mol_scan_dict = generate_scans(
        generate_molecule_isotopologue_lib(peak_props, [1, 2, 3], trivial_names),
        peak_props,
        generate_interval_tree(peak_props),
        NucleosideFragmentor(),
        UniformNoiseInjector(),
        {"gradient_length": 0.06, "min_intensity": 0, "isolation_window_width": 0.2},
    )
    assert len(scans) == 1  # one MS1 with related MS2 scans
    assert len(scans[0]) == 2  # first element is a MS1 scan and a list of MS2 scans
    assert len(scans[0][1]) == 1  # list of MS2 scans has one scan


def test_generate_molecule_isotopologue_lib():
    peak_props = {
        "uridine": {
            "charge": 2,
            "chemical_formula": "+C(9)H(11)N(2)O(6)",
            "trivial_name": "uridine",
            "scan_start_time": 0,
            "peak_width": 0.06,  # peak as long as ms_rt_diff and gradient length
            "peak_function": "gauss",
            "peak_params": {"sigma": 1},  # 10% of peak width,
            "peak_scaling_factor": 1,
        }
    }
    charges = [1, 2, 3]
    trivial_names = {"+C(9)H(11)N(2)O(6)": "uridine"}
    lib = generate_molecule_isotopologue_lib(peak_props, charges, trivial_names)
    assert list(lib.keys()) == ["uridine"]
    assert np.isclose(
        lib["uridine"]["mz"],
        [
            122.53813200786999,
            123.03962222482085,
            123.54123939611428,
            124.04272136528171,
            124.5442376065953,
            125.04387103562883,
            125.54630817011011,
        ],
    ).all()
    assert np.isclose(
        lib["uridine"]["i"],
        [
            1.0,
            0.10819885003176334,
            0.01762911869579814,
            0.0013546885053562444,
            6.246133599652487e-05,
            4.563511253647458e-07,
            4.800737885333426e-10,
        ],
    ).all()


def test_write_mzml_one_spec():
    tempfile = NamedTemporaryFile("wb")
    peak_props = {
        "uridine": {
            "charge": 2,
            "chemical_formula": "+C(9)H(11)N(2)O(6)",
            "trivial_name": "uridine",
            "scan_start_time": 0,
            "peak_width": 0.03,  # peak as long as ms_rt_diff and gradient length
            "peak_function": "gauss",
            "peak_params": {"sigma": 1},  # 10% of peak width,
            "peak_scaling_factor": 1,
        }
    }
    fragmentor = NucleosideFragmentor()
    noise = UniformNoiseInjector()
    mzml_params = {"gradient_length": 0.03, "ms_rt_diff": 0.03}

    write_mzml(tempfile, peak_props, fragmentor, noise, mzml_params)
    reader = pymzml.run.Reader(tempfile.name)
    assert reader.get_spectrum_count() == 1


def test_dynacmic_exclusion():
    tempfile = NamedTemporaryFile("wb")
    peak_props = {
        "uridine": {
            "charge": 2,
            "chemical_formula": "+C(9)H(11)N(2)O(6)",
            "trivial_name": "uridine",
            "scan_start_time": 0,
            "peak_width": 5,  # peak as long as ms_rt_diff and gradient length
            "peak_function": "gauss",
            "peak_params": {"sigma": 1},  # 10% of peak width,
            "peak_scaling_factor": 1e5,
        }
    }
    trivial_names = {"+C(9)H(11)N(2)O(6)": "uridine"}
    # dynamic_exclusion is bigger than gradient length, so we expect only one MS2 fragment spectrum
    default_mzml_params = {
        "gradient_length": 5,
        "min_intensity": 100,
        "isolation_window_width": 0.5,
        "ion_target": 3e6,
        "ms_rt_diff": 0.03,
        "dynamic_exclusion": 30,  # in seconds
    }
    scans, mol_scan_dict = generate_scans(
        generate_molecule_isotopologue_lib(peak_props, [1, 2, 3], trivial_names),
        peak_props,
        generate_interval_tree(peak_props),
        NucleosideFragmentor(),
        UniformNoiseInjector(),
        default_mzml_params,
    )
    number_fragment_specs = 0
    for ms1, ms2_list in scans:
        for ms2 in ms2_list:
            if len(ms2["mz"]) > 0:
                number_fragment_specs += 1
    # fails when dynamic_exclusion is set to 0
    # since there will be 15 fragments specs
    assert number_fragment_specs == 1

    # breakpoint()
