"""Callables for fragmenting molecules.

Upon calling the callabe, a list/np.array of mz and intensities should be returned.
Arguments should be passed via *args and **kwargs
"""
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Union
import sys
import shutil
import subprocess
import os
import csv

import numpy as np
import pandas as pd
import pyqms
from loguru import logger
from pyteomics import mass

import smiter
from peptide_fragmentor import PeptideFragment0r
from smiter.ext.nucleoside_fragment_kb import (
    KB_FRAGMENTATION_INFO as pyrnams_nucleoside_fragment_kb,
)
from smiter.lib import calc_mz

try:
    from smiter.ext.nucleoside_fragment_kb import KB_FRAGMENTATION_INFO
except ImportError:  # pragma: no cover
    print("Nucleoside fragmentation KB not available")  # pragma: no cover


class AbstractFragmentor(ABC):
    """Summary."""

    @abstractmethod
    def __init__(self):
        """Summary."""
        pass  # pragma: no cover

    @abstractmethod
    def fragment(self, entity):
        """Summary.

        Args:
            entity (TYPE): Description
        """
        pass  # pragma: no cover


class PeptideFragmentor(AbstractFragmentor):
    """Summary."""

    def __init__(self, *args, **kwargs):
        """Summary."""
        logger.info("Initialize PeptideFragmentor")
        self.args = args
        self.kwargs = kwargs
        self.fragger = PeptideFragment0r()

    # @profile
    def fragment(self, entities):
        """Summary.

        Args:
            entity (TYPE): Description
        """
        if isinstance(entities, str):
            entities = [entities]
        frames = []
        for entity in entities:
            # logger.debug(f"Fragment {entity}")
            results_table = self.fragger.fragment(entity, **self.kwargs)
            frames.append(results_table)
        final_table = pd.concat(frames)
        i = np.array([100 for i in range(len(final_table))])
        mz_i = np.stack((final_table["mz"], i), axis=1)
        return mz_i


class PeptideFragmentorPyteomics(AbstractFragmentor):
    def __init__(self, *args, **kwargs):
        pass

    # @profile
    def _fragments(self, peptide, types=("b", "y"), maxcharge=1):
        for i in range(1, len(peptide) - 1):
            for ion_type in types:
                for charge in range(1, maxcharge + 1):
                    if ion_type[0] in "abc":
                        yield mass.fast_mass(
                            peptide[:i], ion_type=ion_type, charge=charge
                        )
                    else:
                        yield mass.fast_mass(
                            peptide[i:], ion_type=ion_type, charge=charge
                        )

    # @profile
    def fragment(self, entities):
        mz = []
        for e in entities:
            a = list(self._fragments(e))
            mz.extend(a)
        i = np.array([100 for x in range(len(mz))])
        mz = np.array(mz)
        mz_i = np.stack((mz, i), axis=1)
        return mz_i


class NucleosideFragmentor(AbstractFragmentor):
    """Summary."""

    def __init__(
        self,
        nucleotide_fragment_kb: Dict[str, dict] = None,
        raise_error_for_non_existing_fragments=True,
    ):
        """Summary."""
        logger.info("Initialize NucleosideFragmentor")
        if nucleotide_fragment_kb is None:
            nucleoside_fragment_kb = pyrnams_nucleoside_fragment_kb
        self.raise_error_for_non_existing_fragments = (
            raise_error_for_non_existing_fragments
        )
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

    def fragment(
        self, entities: Union[list, str], raise_error_for_non_existing_fragments=False
    ):
        """Summary.

        Args:
            entity (TYPE): Description
        """
        if isinstance(entities, str):
            entities = [entities]
        m = []
        for entity in entities:
            if raise_error_for_non_existing_fragments is True:
                masses = self.nuc_to_fragments[entity]
            else:
                masses = self.nuc_to_fragments.get(entity, [])
            m.extend(masses)
            # logger.debug(masses)
            # should overlapping peaks be divided into two very similar ones?
        m = sorted(list(set(m)))
        # logger.debug(m)
        return np.array([(mass, 1) for mass in m])


class LipidFragmentor(AbstractFragmentor):
    """Summary."""

    def __init__(
        self,
        lipid_input_csv: str = None,
        raise_error_for_non_existing_fragments=True,
    ):
        """Use LipidCreator to calculate precursor transitions of lipids."""
        self.lip_to_fragments = {}
        # TODO run lipid fragmenter here, read output file and collect results in dict
        commands: List[str] = []
        if sys.platform == "linux" or sys.platform == "darwin":
            commands.append("mono")
            lipid_creator_path = shutil.which("LipidCreator.exe")
        else:
            # will this work under windows?
            lipid_creator_path = "LipidCreator"
        commands.extend(
            [lipid_creator_path, "transitionlist", lipid_input_csv, "lipid_output.csv"]
        )
        proc = subprocess.run(commands)
        with open("lipid_output.csv") as fin:
            for line in csv.DictReader(fin):
                if line["PrecursorName"] not in self.lip_to_fragments:
                    self.lip_to_fragments[line["PrecursorName"]] = []
                self.lip_to_fragments[line["PrecursorName"]].append(
                    float(line["ProductMz"])
                )
        os.remove("lipid_output.csv")

    def fragment(
        self, entities: Union[list, str], raise_error_for_non_existing_fragments=False
    ):
        """Summary.

        Args:
            entity (TYPE): Description
        """
        if isinstance(entities, str):
            entities = [entities]
        m = []
        for entity in entities:
            if raise_error_for_non_existing_fragments is True:
                masses = self.lip_to_fragments[entity]
            else:
                masses = self.lip_to_fragments.get(entity, [])
            m.extend(masses)
            # logger.debug(masses)
            # should overlapping peaks be divided into two very similar ones?
        m = sorted(list(set(m)))
        # logger.debug(m)
        return np.array([(mass, 1) for mass in m])
