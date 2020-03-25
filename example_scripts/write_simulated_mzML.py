#!/usr/bin/env python3
import numpy as np
import pymzml
from pymzml.plot import Factory

from smiter.fragmentation_functions import NucleosideFragmentor
from smiter.synthetic_mzml import write_mzml


def main():
    molecules = ["+C(10)H(12)N(4)O(5)"]
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "charge": 1,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},
            "trivial_name": "inosine",
        }
    }
    mzml_params = {"gradient_length": 30, "ms_rt_diff": 0.01}
    file = "inosine.mzML"
    fragmentor = NucleosideFragmentor()
    mzml_path = write_mzml(file, molecules, fragmentor, peak_props, mzml_params)


if __name__ == "__main__":
    main()
