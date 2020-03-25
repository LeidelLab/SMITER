#!/usr/bin/env python3

from tempfile import NamedTemporaryFile

import pymzml
from pymzml.plot import Factory
from smiter.synthetic_mzml import write_mzml
from smiter.fragmentation_functions import NucleosideFragmentor


def main():
    molecules = ["+C(10)H(12)N(4)O(5)"]
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "trivial_name": 'inosine',
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 0.1},
        }
    }
    mzml_params = {
        "gradient_length": 30,
    }

    fragmentor = NucleosideFragmentor()

    file = NamedTemporaryFile("wb")
    mzml_path = write_mzml(file, molecules, fragmentor, peak_props, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    tic_tuples = []
    for pos, spec in enumerate(reader):
        if spec.ms_level == 1:
            # print(spec.TIC)
            tic_tuples.append((pos, spec.TIC))

    pf = Factory()
    pf.new_plot()
    pf.add(tic_tuples, color=(0, 0, 0), style="lines", title="TIC")
    pf.save(
        "example_chromatogram.html",
        layout={"xaxis": {"title": "Retention time"}, "yaxis": {"title": "TIC"}},
    )


if __name__ == "__main__":
    main()
