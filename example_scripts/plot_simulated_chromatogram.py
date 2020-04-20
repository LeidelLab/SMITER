#!/usr/bin/env python3

from tempfile import NamedTemporaryFile

import pymzml
from pymzml.plot import Factory

from smiter.fragmentation_functions import NucleosideFragmentor
from smiter.synthetic_mzml import write_mzml


def main():
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "trivial_name": "inosine",
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},
        }
    }
    mzml_params = {
        "gradient_length": 30,
    }

    fragmentor = NucleosideFragmentor()

    file = NamedTemporaryFile("wb")
    mzml_path = write_mzml(file,  peak_props, fragmentor, mzml_params)
    reader = pymzml.run.Reader(mzml_path)
    tic_tuples = []
    for pos, spec in enumerate(reader):
        if spec.ms_level == 1:
            scan_time, scan_unit = spec.scan_time
            if scan_unit == "second":
                scan_time /= 60
            tic_tuples.append((scan_time, spec.TIC))

    pf = Factory()
    pf.new_plot()
    pf.add(tic_tuples, color=(0, 0, 0), style="lines", title="TIC")
    pf.save(
        "example_chromatogram.html",
        layout={
            "xaxis": {"title": "Retention time [minutes]"},
            "yaxis": {"title": "TIC"},
        },
    )


if __name__ == "__main__":
    main()
