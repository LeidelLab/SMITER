"""Core functionality."""
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
