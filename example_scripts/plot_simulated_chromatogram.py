#!/usr/bin/env python3

from tempfile import NamedTemporaryFile

import pymzml
from pymzml.plot import Factory

from smiter.fragmentation_functions import NucleosideFragmentor
from smiter.noise_functions import UniformNoiseInjector, JamssNoiseInjector
from smiter.synthetic_mzml import write_mzml


def main():
    peak_props = {
        "inosine": {
            "trivial_name": "inosine",
            "chemical_formula": "+C(10)H(12)N(4)O(5)",
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss_tail",
            "peak_params": {"sigma": 3},
            "peak_scaling_factor": 1e6,
        },
        "+C(9)H(12)N(2)O(6)": {
            "trivial_name": "pseudouridine",
            "charge": 2,
            "scan_start_time": 30,
            "peak_width": 30,  # seconds
            "peak_function": "gauss_tail",
            "peak_scaling_factor": 2e6,
            "peak_params": {"sigma": 3},
        },
        "+C(10)H(13)N(5)O(7)": {
            "trivial_name": "spiroiminodihydantoin",
            "charge": 2,
            "scan_start_time": 60,
            "peak_width": 30,  # seconds
            "peak_function": "gauss_tail",
            "peak_params": {"sigma": 3},
            "peak_scaling_factor": 3e6,
        }

    }
    mzml_params = {
        "gradient_length": 30,
        "min_intensity": 1,
    }

    fragmentor = NucleosideFragmentor()
    noise_injector = UniformNoiseInjector()
    noise_injector = UniformNoiseInjector(
        dropout=0.0, ppm_noise=0, intensity_noise=0
    )
    noise_injector = JamssNoiseInjector()


    file = NamedTemporaryFile("wb")
    mzml_path = write_mzml(file, peak_props, fragmentor, noise_injector, mzml_params)
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
