import csv
from pprint import pprint

import click

import smiter
import warnings

warnings.simplefilter("ignore")


@click.command()
@click.argument("input_csv")
@click.argument("output_mzml")
def main(input_csv, output_mzml):
    peak_properties = smiter.lib.csv_to_peak_properties(input_csv)
    fragmentor = smiter.fragmentation_functions.NucleosideFragmentor()
    noise_injector = smiter.noise_functions.UniformNoiseInjector(
        dropout=0.3, ppm_noise=5e-6, intensity_noise=0.2
    )
    mzml_params = {"gradient_length": 60}
    smiter.synthetic_mzml.write_mzml(
        output_mzml, peak_properties, fragmentor, noise_injector, mzml_params
    )



if __name__ == "__main__":
    main()
