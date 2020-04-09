"""Core functionality."""
from smiter.params.default_params import default_mzml_params, default_peak_properties

PROTON = 1.00727646677


def calc_mz(mass: float, charge: int):
    """Calculate m/z.

    Args:
        mass (TYPE): Description
        charge (TYPE): Description
    """
    mass = float(mass)
    charge = int(charge)
    calc_mz = (mass + (charge * PROTON)) / charge
    return calc_mz


def check_mzml_params(mzml_params: dict) -> dict:
    """Summary.

    Args:
        mzml_params (dict): Description

    Returns:
        dict: Description

    Raises:
        Exception: Description
    """
    for default_param, default_value in default_mzml_params.items():
        # param not set and default param required
        if (mzml_params.get(default_param, None) is None) and (default_value is None):
            raise Exception(f"mzml parameter {default_param} is required by not set!")
        elif mzml_params.get(default_param, None) is None:
            mzml_params[default_param] = default_value
    return mzml_params


def check_peak_properties(peak_properties: dict) -> dict:
    """Summary.

    Args:
        peak_properties (dict): Description

    Returns:
        dict: Description

    Raises:
        Exception: Description
    """
    for mol, properties in peak_properties.items():
        for default_param, default_value in default_peak_properties.items():
            if (properties.get(default_param, None) is None) and (
                default_value is None
            ):
                raise Exception(
                    f"mzml parameter {default_param} is required by not set!"
                )
            elif properties.get(default_param, None) is None:
                properties[default_param] = default_value
    return peak_properties
