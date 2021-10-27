import csv
from pprint import pprint

import click
import matplotlib.pyplot as plt
import pymzml

import smiter
import warnings

warnings.simplefilter("ignore")


@click.command()
@click.argument("input_csv")
@click.argument("output_mzml")
def main(input_csv, output_mzml):
    peak_properties = smiter.lib.csv_to_peak_properties(input_csv)
    pprint(peak_properties)
    fragmentor = smiter.fragmentation_functions.NucleosideFragmentor()
    noise_injector = smiter.noise_functions.UniformNoiseInjector(
        dropout=0.3, ppm_noise=5e-6, intensity_noise=0.05
    )
    # noise_injector = smiter.noise_functions.UniformNoiseInjector(
    #     dropout=0.0, ppm_noise=0, intensity_noise=0
    # )

    mzml_params = {"gradient_length": 60}
    smiter.synthetic_mzml.write_mzml(
        output_mzml, peak_properties, fragmentor, noise_injector, mzml_params
    )

    rt = []
    i  = []
    with pymzml.run.Reader(output_mzml) as run:
        for spec in run:
            if spec.ms_level == 1:
                rt.append(spec.scan_time_in_minutes() * 60)
                i.append(spec.i.sum())
    plt.plot(rt, i)
    plt.show()

if __name__ == "__main__":
    main()
