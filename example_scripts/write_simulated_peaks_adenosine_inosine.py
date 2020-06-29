#!/usr/bin/env python3
"""Write file containing inosine and adenosine peak."""
import os

import click
import matplotlib.pyplot as plt
import numpy as np
import pymzml
import warnings

from smiter.fragmentation_functions import AbstractFragmentor, NucleosideFragmentor
from smiter.noise_functions import UniformNoiseInjector, JamssNoiseInjector
from smiter.synthetic_mzml import write_mzml

warnings.simplefilter("ignore")


@click.command()
@click.argument("file_path")
def main(file_path):
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
    fragmentor = NucleosideFragmentor()
    # noise_injector= UniformNoiseInjector()
    noise_injector = JamssNoiseInjector()
    with open(file_path, "wb") as fin:
        mzml_path = write_mzml(fin, peak_props, fragmentor, noise_injector, mzml_params)

    reader = pymzml.run.Reader(file_path)

    ino_rt = []
    ino_intensities = []
    adeno_rt = []
    adeno_intensities = []

    ino_mono = 244.06898754897
    adeno_mono = 244.06898754897

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
    plt.plot(ino_rt, ino_intensities)
    plt.plot(adeno_rt, adeno_intensities)
    plt.savefig(os.path.splitext(file_path)[0])


if __name__ == "__main__":
    main()
