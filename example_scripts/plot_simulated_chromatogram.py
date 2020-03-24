#!/usr/bin/env python3

import pymzml
from pymzml.plot import Factory

from tempfile import NamedTemporaryFile

from smiter.synthetic_mzml import write_mzml


def simple_fragmentation_function(mol):
    """Test frag function.

    Args:
        mol (str): molecule to fragment

    Returns:
        list: list of fragment mzs
    """
    return [20]

def main():
    molecules = [
        "+C(10)H(12)N(4)O(5)"
    ]
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {
                "sigma": 0.1
            }
        }
    }
    file = NamedTemporaryFile("wb")
    mzml_path = write_mzml(
        file,
        molecules,
        simple_fragmentation_function,
        peak_props
    )
    reader = pymzml.run.Reader(mzml_path)
    tic_tuples = []
    for pos,spec in enumerate(reader):
        if spec.ms_level == 1:
            # print(spec.TIC)
            tic_tuples.append((pos, spec.TIC))

    pf = Factory()
    pf.new_plot()
    pf.add(
        tic_tuples,
        color=(0, 0, 0),
        style="lines",
        title='TIC'
    )
    pf.save(
        "example_chromatogram.html",
        layout={
            "xaxis": {
                "title": "Retention time"
            },
            "yaxis": {
                "title": "TIC"
            }
        },
    )

if __name__ == "__main__":
    main()