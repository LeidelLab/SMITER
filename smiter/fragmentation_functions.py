"""Callables for fragmenting molecules.

Upon calling the callabe, a list/np.array of mz and intensities should be returned.
Arguments should be passed via *args and **kwargs
"""
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple

import numpy as np
import pyqms

from smiter.ext.nucleoside_fragment_kb import (
    KB_FRAGMENTATION_INFO as pyrnams_nucleoside_fragment_kb,
)
from smiter.lib import calc_mz

try:
    from smiter.ext.nucleoside_fragment_kb import KB_FRAGMENTATION_INFO
except ImportError:
    print("Nucleoside fragmentation KB not available")


class AbstractFragmentor(ABC):
    """Summary."""

    def __init__(self):
        """Summary."""
        super().__init__()

    @abstractmethod
    def fragment(self, entity):
        """Summary.

        Args:
            entity (TYPE): Description
        """
        pass


class PeptideFragmentor(AbstractFragmentor):
    """Summary."""

    def __init__(self):
        """Summary."""
        pass

    def fragment(self, entity):
        """Summary.

        Args:
            entity (TYPE): Description
        """
        pass


class NucleosideFragmentor(AbstractFragmentor):
    """Summary."""

    def __init__(self, nucleotide_fragment_kb: Dict[str, dict] = None):
        """Summary."""
        if nucleotide_fragment_kb is None:
            nucleoside_fragment_kb = pyrnams_nucleoside_fragment_kb

        nuc_to_fragments: Dict[str, List[float]] = {}
        cc = pyqms.chemical_composition.ChemicalComposition()
        for nuc_name, nuc_dict in nucleoside_fragment_kb.items():
            nuc_to_fragments[nuc_name] = []
            for frag_name, frag_cc_dict in nucleoside_fragment_kb[nuc_name][
                "fragments"
            ].items():
                cc.use(f"+{frag_cc_dict['formula']}")
                m = cc._mass()
                nuc_to_fragments[nuc_name].append(calc_mz(m, 1))
        self.nuc_to_fragments = nuc_to_fragments

    def fragment(self, entity):
        """Summary.

        Args:
            entity (TYPE): Description
        """
        masses = self.nuc_to_fragments[entity]
        return np.array([(mass, 100) for mass in masses])
