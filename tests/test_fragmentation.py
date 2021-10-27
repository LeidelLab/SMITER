"""Summary."""
import os
import numpy as np
import pytest

import smiter
from smiter.fragmentation_functions import (
    NucleosideFragmentor,
    PeptideFragmentor,
    LipidFragmentor,
)


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
    expected_mzs = np.array([119.03522254717, 136.0617716478])
    assert ((peaks[:, 0] - expected_mzs) < 0.001).all()


@pytest.mark.skipif("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true", reason="Skipping this test on Travis CI.")
def test_fragment_lipid():
    test_lipid_file = "test_lipids.txt"
    test_lipid_file = os.path.abspath(test_lipid_file)
    with open(test_lipid_file, "wt") as fout:
        fout.write("PC 18:0/12:0\n")
        fout.write("PE 18:3;1-16:2")
    fragger = LipidFragmentor(test_lipid_file)
    masses = fragger.fragment("PC 18:0/12:0")
    assert np.allclose(masses, np.array([184.0733]))
