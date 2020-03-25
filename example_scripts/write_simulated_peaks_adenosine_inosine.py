#!/usr/bin/env python3
"""Write file containing inosine and adenosine peak."""
import click
import numpy as np

from smiter.fragmentation_functions import AbstractFragmentor, NucleosideFragmentor
from smiter.synthetic_mzml import write_mzml


@click.command()
@click.argument("file_path")
def main(file_path):
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
    fragmentor = NucleosideFragmentor()
    with open(file_path, "wb") as fin:
        mzml_path = write_mzml(fin, molecules, fragmentor, peak_props, mzml_params)


if __name__ == "__main__":
    main()
