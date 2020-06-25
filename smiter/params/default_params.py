# None means required, other is default value
# non existent params are optional
default_peak_properties = {
    # trivial_name
    "charge": 2,
    "scan_start_time": None,  # required
    "peak_width": None,  # required
    # "peak_function": "gauss",
}

default_mzml_params = {"gradient_length": None, "isolation_window_width": 0.5}
