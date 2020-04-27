#!/usr/bin/env python3
"""Write file containing inosine and adenosine peak."""
import os

import click
import matplotlib.pyplot as plt
import numpy as np
import pymzml

from smiter.fragmentation_functions import AbstractFragmentor, NucleosideFragmentor
from smiter.noise_functions import UniformNoiseInjector
from smiter.synthetic_mzml import write_mzml


@click.command()
@click.argument("file_path")
def main(file_path):
    peak_props = {
        "+C(10)H(12)N(4)O(5)": {
            "charge": 2,
            "trivial_name": "inosine",
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
            "peak_scaling_factor": 0.5 * 1e6,
        },
        "+C(10)H(13)N(5)O(4)": {
            "trivial_name": "adenosine",
            "charge": 2,
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
    fragmentor = NucleosideFragmentor()
    noise_injector= UniformNoiseInjector()
    with open(file_path, "wb") as fin:
        mzml_path = write_mzml(fin, peak_props, fragmentor, noise_injector, mzml_params)

    reader = pymzml.run.Reader(file_path)

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
    plt.plot(ino_rt, ino_intensities)
    plt.plot(adeno_rt, adeno_intensities)
    plt.savefig(os.path.splitext(file_path)[0])


if __name__ == "__main__":
    main()
