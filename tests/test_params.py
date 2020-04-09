"""Summary.
"""
import pytest

from smiter.lib import check_mzml_params, check_peak_properties


def test_check_mzml_params():
    mzml_params = {
        "gradient_length": 45,
    }
    mzml_params = check_mzml_params(mzml_params)
    assert mzml_params["gradient_length"] == 45


def test_check_faulty_mzml_params():
    mzml_params = {}
    with pytest.raises(Exception):
        mzml_params = check_mzml_params(mzml_params)


def test_check_peak_properties_missing_default():
    peak_props = {
        "ELVISLIVES": {
            "trivial_name": "ELVISLIVES",
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width
        },
    }
    peak_props = check_peak_properties(peak_props)
    assert peak_props["ELVISLIVES"]["charge"] == 2


def test_check_peak_properties_missing_opt():
    peak_props = {
        "ELVISLIVES": {
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width
        },
    }
    peak_props = check_peak_properties(peak_props)
    # sets charge but leaves out trivial_name
    assert peak_props["ELVISLIVES"]["charge"] == 2
    assert "trivial_name" not in peak_props["ELVISLIVES"]


def test_check_peak_properties_missing_required():
    peak_props = {
        "ELVISLIVES": {
            "trivial_name": "ELVISLIVES",
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width
        },
    }
    with pytest.raises(Exception):
        peak_props = check_peak_properties(peak_props)
