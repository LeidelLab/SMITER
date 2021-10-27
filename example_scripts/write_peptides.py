#!//usr/bin/env python
from tempfile import NamedTemporaryFile
import smiter
from smiter.fragmentation_functions import (
    AbstractFragmentor,
    NucleosideFragmentor,
    PeptideFragmentor,
)
from smiter.noise_functions import GaussNoiseInjector, UniformNoiseInjector
from smiter.synthetic_mzml import write_mzml


def main():
    file = NamedTemporaryFile("wb")
    peak_props = {
        "ELVISLIVES": {
            "trivial_name": "ELVISLIVES",
            "charge": 2,
            "scan_start_time": 0,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width
        },
        "ELVISLIVSE": {
            "charge": 2,
            "trivial_name": "ELVISLIVSE",
            "scan_start_time": 15,
            "peak_width": 30,  # seconds
            "peak_function": "gauss",
            "peak_params": {"sigma": 3},  # 10% of peak width,
        },
    }
    mzml_params = {
        "gradient_length": 45,
    }
    fragmentor = PeptideFragmentor()
    noise_injector = GaussNoiseInjector(variance=0.05)
    mzml_path = write_mzml(file, peak_props, fragmentor, noise_injector, mzml_params)


if __name__ == '__main__':
    main()
