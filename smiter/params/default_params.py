# None means required, other is default value
# non existent params are optional
default_peak_properties = {
    # trivial_name
    "charge": 2,
    "scan_start_time": None,  # required
    "peak_width": None,  # required
    # "peak_function": "gauss",
}

default_mzml_params = {
    "gradient_length": None,
    "min_intensity": 100,
    "max_intensity": 1e10,  # what would be reasonable?
    "isolation_window_width": 0.5,
    "ion_target": 3e6,
    "ms_rt_diff": 0.03,
    "dynamic_exclusion": 30,  # in seconds
    "max_ms2_spectra": 10,
    "mz_lower_limit": 100,
    "mz_upper_limit": 1600,
}
