"""Summary."""
import numpy as np

import smiter
from smiter.fragmentation_functions import NucleosideFragmentor, PeptideFragmentor


def test_fragment_peptide():
    """Summary."""
    fragger = PeptideFragmentor(charges=1, ions=["y"])
    peaks = fragger.fragment("L")
    expected_mzs = np.array([132.101])
    assert ((peaks[:, 0] - expected_mzs) < 0.001).all()


def test_fragment_nucleotide():
    """Summary."""
    fragger = NucleosideFragmentor()
    peaks = fragger.fragment("adenosine")
    expected_mzs = np.array([136.0617716478, 119.03522254717])
    assert ((peaks[:, 0] - expected_mzs) < 0.001).all()
