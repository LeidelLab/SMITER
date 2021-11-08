=====
Usage
=====

To use SMITER in a project::

    1. create the peak properties dict (You can convert a csv file (as in example_data) to csv using `smiter.lib.csv_to_peak_properties`)
    2. Choose a fragmentor (e.g. fragmentation_functions.PeptideFragmentor)
    3. Choose a Noise generator (e.g. smiter.noise_functions.UniformNoiseInjector)
    4. Define general params (e.g. gradient_length)
    5. Run the simulation and write the resulting mzML using `smiter.synthetic_mzml.write_mzml` with the peak_properties, noise injector, fragmentor etc. as parameters, see the API documentation.


When in doubt, check the complex_mix_simulation.py script