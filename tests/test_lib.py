"""Summary.
"""
import csv
import os
from tempfile import NamedTemporaryFile

import pytest

from smiter.lib import (
    check_mzml_params,
    check_peak_properties,
    csv_to_peak_properties,
    peak_properties_to_csv,
)


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


def test_csv_to_peak_properties():
    csv_file = os.path.join(os.path.dirname(__file__), "data", "molecules_test.csv")
    peak_properties = csv_to_peak_properties(csv_file)
    assert "2′-O-methylcytidine" in peak_properties
    assert (
        peak_properties["2′-O-methylcytidine"]["trivial_name"] == "2′-O-methylcytidine"
    )
    assert peak_properties["2′-O-methylcytidine"]["scan_start_time"] == 10
    assert peak_properties["2′-O-methylcytidine"]["peak_scaling_factor"] == 1e6
    # default
    assert peak_properties["2′-O-methylcytidine"]["charge"] == 2

    assert "5-methoxycarbonylmethyluridine" in peak_properties
    assert (
        peak_properties["5-methoxycarbonylmethyluridine"]["trivial_name"]
        == "5-methoxycarbonylmethyluridine"
    )
    assert peak_properties["5-methoxycarbonylmethyluridine"]["scan_start_time"] == 12
    assert (
        peak_properties["5-methoxycarbonylmethyluridine"]["peak_scaling_factor"] == 2e6
    )
    # default
    assert peak_properties["5-methoxycarbonylmethyluridine"]["charge"] == 2


def test_peak_properties_to_csv():
    peak_properties = {
        "2′-O-methylcytidine": {
            "trivial_name": "2′-O-methylcytidine",
            "chemical_formula": "+C(10)H(15)N(3)O(5)",
            "charge": 2,
            "scan_start_time": 10.0,
            "peak_function": "gauss",
            "peak_params": {"sigma": 2},
            "peak_scaling_factor": 1000000.0,
            "peak_width": 30,
        },
        "5-methoxycarbonylmethyluridine": {
            "trivial_name": "5-methoxycarbonylmethyluridine",
            "chemical_formula": "+C(12)H(16)N(2)O(8)",
            "charge": 2,
            "scan_start_time": 12.0,
            "peak_function": "gauss",
            "peak_params": {"sigma": 2},
            "peak_scaling_factor": 2000000.0,
            "peak_width": 30,
        },
    }
    csv_file = "out.csv"
    fname = peak_properties_to_csv(peak_properties, csv_file)
    assert fname == csv_file
    with open(csv_file) as fout:
        reader = csv.DictReader(fout)
        lines = [l for l in reader]
    print(lines)

    assert lines[0]["chemical_formula"] == "+C(10)H(15)N(3)O(5)"
    assert lines[0]["scan_start_time"] == "10.0"
    assert lines[0]["peak_function"] == "gauss"
    assert lines[0]["peak_params"] == "sigma=2"
    assert lines[0]["peak_scaling_factor"] == "1000000.0"
    assert lines[0]["peak_width"] == "30"
    # default
    assert lines[0]["charge"] == "2"

    assert lines[1]["chemical_formula"] == "+C(12)H(16)N(2)O(8)"
    assert lines[1]["scan_start_time"] == "12.0"
    assert lines[1]["peak_function"] == "gauss"
    assert lines[1]["peak_params"] == "sigma=2"
    assert lines[1]["peak_scaling_factor"] == "2000000.0"
    assert lines[1]["peak_width"] == "30"
    # default
    assert lines[0]["charge"] == "2"
