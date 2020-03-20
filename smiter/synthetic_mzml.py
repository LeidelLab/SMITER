"""Main module."""
import io
from typing import Callable, Dict, List, Tuple, Union

import numpy as np

import pyqms
import scipy as sci
from psims.mzml import MzMLWriter
from smiter.peak_distribution import distributions


class Scan(dict):
    """Summary."""

    def __init__(self, data: dict = None):
        """Summary.

        Args:
            dict (TYPE): Description
        """
        if data is not None:
            self.update(data)

    @property
    def mz(self):
        """Summary."""
        v = self.get("mz", None)
        return v

    @property
    def i(self):
        """Summary."""
        v = self.get("i", None)
        return v

    @i.setter
    def i(self, i):
        """Summary."""
        self["i"] = i

    @property
    def id(self):
        """Summary."""
        v = self.get("id", None)
        return v

    @property
    def precursor_mz(self):
        """Summary."""
        v = self.get("precursor_mz", None)
        return v

    @property
    def precursor_i(self):
        """Summary."""
        v = self.get("precursor_i", None)
        return v

    @property
    def precursor_charge(self):
        """Summary."""
        v = self.get("precursor_charge", None)
        return v

    @property
    def precursor_scan_id(self):
        """Summary."""
        v = self.get("precursor_scan_id", None)
        return v

    @property
    def retention_time(self):
        """Summary."""
        v = self.get("rt", None)
        return v

    @property
    def ms_level(self):
        v = self.get("ms_level", None)
        return v


def write_mzml(
    file: Union[str, io.TextIOWrapper],
    molecules: List[str],
    fragmentation_function: Callable[[str], List[Tuple[float, float]]],
    peak_properties: Dict[str, dict] = None,
) -> str:
    """Write mzML file with chromatographic peaks and fragment spectra for the given molecules.

    Args:
        file (Union[str, io.TextIOWrapper]): Description
        molecules (List[str]): Description
        fragmentation_function (Callable[[str], List[Tuple[float, float]]], optional): Description
        peak_properties (Dict[str, dict], optional): Description
    """
    filename = file if isinstance(file, str) else file.name
    scans = []
    isotopologue_lib = generate_molecule_isotopologue_lib(molecules)
    if len(isotopologue_lib) > 0:
        scans, scan_dict = generate_scans(
            isotopologue_lib, peak_properties, fragmentation_function
        )
    else:
        scans, scan_dict = [], {}
    # TODO also rescale ms2 scans?
    # TODO ms2 scan intensitity as fraction of precursor intensity?
    scans = rescale_ms1_scans(scans, scan_dict, peak_properties=peak_properties)
    write_scans(file, scans)
    return filename


def rescale_ms1_scans(scans: list, scan_dict: dict, peak_properties=None):
    """Summary

    Args:
        scans (list): Description
        scan_dict (dict): Description
    """
    if peak_properties is None:
        peak_properties = {}
    for molecule in scan_dict:
        scan_num = len(scan_dict[molecule]["ms1_scans"])
        scale_func = peak_properties[f"+{molecule}"]["peak_function"]
        for i, scan_list in enumerate(scans):
            s = scan_list[0]
            if scale_func == "gauss":
                scale_factor = distributions[scale_func](
                    i,
                    mu=scan_num / 2,
                    sigma=scan_num
                    * peak_properties[f"+{molecule}"]["peak_params"]["sigma"],
                )
            elif scale_func == "gamma":
                scale_factor = distributions[scale_func](
                    i,
                    a=peak_properties[f"+{molecule}"]["peak_params"]["a"],
                    scale=peak_properties[f"+{molecule}"]["peak_params"]["scale"],
                )
            elif scale_func is None:
                scale_factor = 1
            s.i *= scale_factor * 1e6  # TODO make this magic number a parameter
    return scans


def generate_scans(isotopologue_lib, peak_properties, fragmentation_function):
    gradient_length = 30
    ms_rt_diff = 0.03

    t = 0
    ms1_scans = 0
    ms2_scans = 0

    mol_scan_dict = {}
    scans = []
    i = 0

    if len(isotopologue_lib) == 0:
        return []

    mol_scan_dict = {
        mol: {"ms1_scans": [], "ms2_scans": []} for mol in isotopologue_lib
    }

    while t < gradient_length:
        for mol in isotopologue_lib:
            mol_plus = f"+{mol}"
            if (peak_properties[mol_plus]["scan_start_time"] <= t) and (
                (
                    peak_properties[mol_plus]["scan_start_time"]
                    + peak_properties[mol_plus]["peak_width"]
                )
                >= t
            ):
                mol_scan_dict[mol]["ms1_scans"].append(i)
                s = Scan(
                    {
                        "mz": isotopologue_lib[mol]["mz"],
                        "rt": t,
                        "i": np.array(isotopologue_lib[mol]["i"]),
                        "id": i,
                    }
                )
                prec_scan_id = i
                i += 1
                scans.append((s, []))
            else:
                scans.append((Scan({"mz": [], "rt": t, "i": [], "id": i}), []))
                prec_scan_id = i
                i += 1

            t += ms_rt_diff

            for ms2 in range(10):
                frag_mz = fragmentation_function(mol)
                ms2_scan = Scan(
                    {
                        "mz": frag_mz,
                        "i": np.array([100 for _ in range(len(frag_mz))]),
                        "rt": t,
                        "id": i,
                        "precursor_mz": 100,
                        "precursor_i": 100,
                        "precursor_charge": 1,
                        "precursor_scan_id": prec_scan_id,
                    }
                )
                mol_scan_dict[mol]["ms2_scans"].append(i)

                scans[-1][1].append(ms2_scan)
                t += ms_rt_diff
                i += 1
    return scans, mol_scan_dict


def generate_molecule_isotopologue_lib(molecules):
    """Summary.

    Args:
        molecules (TYPE): Description
    """
    if len(molecules) > 0:
        lib = pyqms.IsotopologueLibrary(molecules=molecules, charges=[1],)
        reduced_lib = {}
        for mol in lib:
            data = lib[mol]["env"][(("N", "0.000"),)]
            reduced_lib[mol] = {"mz": data[1]["mz"], "i": data["relabun"]}
    else:
        reduced_lib = {}
    return reduced_lib


def write_scans(
    file: Union[str, io.TextIOWrapper], scans: List[Tuple[Scan, List[Scan]]]
) -> None:
    """Generate given scans to mzML file.

    Args:
        file (Union[str, io.TextIOWrapper]): Description
        scans (List[Tuple[Scan, List[Scan]]]): Description

    Returns:
        None: Description
    """
    with MzMLWriter(file) as writer:
        # Add default controlled vocabularies
        writer.controlled_vocabularies()
        # Open the run and spectrum list sections
        with writer.run(id="Simulated Run"):
            spectrum_count = len(scans) + sum([len(products) for _, products in scans])
            with writer.spectrum_list(count=spectrum_count):
                for scan, products in scans:
                    # Write Precursor scan
                    writer.write_spectrum(
                        scan.mz,
                        scan.i,
                        id=scan.id,
                        params=[
                            "MS1 Spectrum",
                            {"ms level": 1},
                            {"total ion current": sum(scan.i)},
                        ],
                    )
                    # Write MSn scans
                    for prod in products:
                        prec_info = {"mz": 100}
                        writer.write_spectrum(
                            prod.mz,
                            prod.i,
                            id=prod.id,
                            params=[
                                "MSn Spectrum",
                                {"ms level": 2},
                                {"total ion current": sum(prod.i)},
                            ],
                            # TOFIX adding precursor information makes psims crash
                            # Include precursor information
                            # precursor_information=prec_info
                            # precursor_information={
                            # "mz": prod.precursor_mz,
                            # "intensity": prod.precursor_i,
                            # "charge": prod.precursor_charge,
                            # "scan_id": prod.precursor_scan_id,
                            # },
                        )
    return
