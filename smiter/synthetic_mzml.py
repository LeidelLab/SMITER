"""Main module."""
import io
from typing import Callable, Dict, List, Tuple, Union

import numpy as np
import pyqms
import scipy as sci
from psims.mzml import MzMLWriter

from smiter.fragmentation_functions import AbstractFragmentor
from smiter.lib import calc_mz, check_mzml_params, check_peak_properties
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

    # @i.setter
    # def i(self, i):
    #     """Summary."""
    #     self["i"] = i

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

    # @property
    # def precursor_scan_id(self):
    #     """Summary."""
    #     v = self.get("precursor_scan_id", None)
    #     return v

    @property
    def retention_time(self):
        """Summary."""
        v = self.get("rt", None)
        return v

    # @property
    # def ms_level(self):
    #     """Summary.

    #     Returns:
    #         TYPE: Description
    #     """
    #     v = self.get("ms_level", None)
    #     return v


def write_mzml(
    file: Union[str, io.TextIOWrapper],
    peak_properties: Dict[str, dict],
    # molecules: List[str],
    fragmentor: AbstractFragmentor,
    mzml_params: Dict[str, Union[int, float, str]],
) -> str:
    """Write mzML file with chromatographic peaks and fragment spectra for the given molecules.

    Args:
        file (Union[str, io.TextIOWrapper]): Description
        molecules (List[str]): Description
        fragmentation_function (Callable[[str], List[Tuple[float, float]]], optional): Description
        peak_properties (Dict[str, dict], optional): Description
    """
    # check params and raise Exception(s) if necessary
    mzml_params = check_mzml_params(mzml_params)
    peak_properties = check_peak_properties(peak_properties)

    filename = file if isinstance(file, str) else file.name
    scans = []
    # pass list of all charge states in peak properties
    trivial_names = {
        key: value["trivial_name"] for key, value in peak_properties.items()
    }
    # dicts are sorted, language specification since python 3.7+
    isotopologue_lib = generate_molecule_isotopologue_lib(
        peak_properties, trivial_names=trivial_names
    )
    # if len(isotopologue_lib) > 0:
    scans, scan_dict = generate_scans(
        isotopologue_lib, peak_properties, fragmentor, mzml_params
    )
    # else:
    # scans, scan_dict = [], {}
    # TODO also rescale ms2 scans?
    # TODO ms2 scan intensitity as fraction of precursor intensity?
    # scans = rescale_ms1_scans(
    #     scans,
    #     scan_dict,
    #     peak_properties=peak_properties,
    #     isotopologue_lib=isotopologue_lib,
    # )
    write_scans(file, scans)
    return filename


def rescale_intensity(
    i: float, rt: float, molecule: str, peak_properties: dict, isotopologue_lib: dict
):
    """Rescale intensity value for a given molecule according to scale factor and distribution function.

    Args:
        i (TYPE): Description
        rt (TYPE): Description
        molecule (TYPE): Description
        peak_properties (TYPE): Description
        isotopologue_lib (TYPE): Description

    Returns:
        TYPE: Description
    """
    scale_func = peak_properties[f"{molecule}"]["peak_function"]
    rt_max = (
        peak_properties[f"{molecule}"]["scan_start_time"]
        + peak_properties[f"{molecule}"]["peak_width"]
    )
    mu = (
        peak_properties[f"{molecule}"]["scan_start_time"]
        + 0.5 * peak_properties[f"{molecule}"]["peak_width"]
    )
    if scale_func == "gauss":
        dist_scale_factor = distributions[scale_func](
            # i,
            rt,
            mu=mu,
            sigma=peak_properties[f"{molecule}"]["peak_params"].get(
                "sigma", peak_properties[f"{molecule}"]["peak_width"] / 10
            ),
        )
    elif scale_func == "gamma":
        dist_scale_factor = distributions[scale_func](
            rt,
            a=peak_properties[f"{molecule}"]["peak_params"]["a"],
            scale=peak_properties[f"{molecule}"]["peak_params"]["scale"],
        )
    elif scale_func is None:
        dist_scale_factor = 1
    i *= dist_scale_factor * peak_properties[f"{molecule}"].get(
        "peak_scaling_factor", 1e3
    )
    return i


def generate_scans(
    isotopologue_lib: dict,
    peak_properties: dict,
    fragmentor: AbstractFragmentor,
    mzml_params: dict,
):
    """Summary.

    Args:
        isotopologue_lib (TYPE): Description
        peak_properties (TYPE): Description
        fragmentation_function (A): Description
        mzml_params (TYPE): Description
    """
    gradient_length = mzml_params["gradient_length"]
    ms_rt_diff = mzml_params.get("ms_rt_diff", 0.03)

    t = 0
    ms1_scans = 0
    ms2_scans = 0

    mol_scan_dict: Dict[str, Dict[str, list]] = {}
    scans: List[Tuple[Scan, List[Scan]]] = []
    i = 0

    # if len(isotopologue_lib) == 0:
    #     return []

    mol_scan_dict = {
        mol: {"ms1_scans": [], "ms2_scans": []} for mol in isotopologue_lib
    }
    molecules = list(isotopologue_lib.keys())
    while t < gradient_length:
        scan_peaks: List[Tuple[float, float]] = []
        mol_i = []
        mol_monoisotopic = {}
        for mol in isotopologue_lib:
            mol_plus = f"{mol}"
            if (peak_properties[mol_plus]["scan_start_time"] <= t) and (
                (
                    peak_properties[mol_plus]["scan_start_time"]
                    + peak_properties[mol_plus]["peak_width"]
                )
                >= t
            ):
                mz = isotopologue_lib[mol]["mz"]
                intensity = np.array(isotopologue_lib[mol]["i"])
                intesity = rescale_intensity(
                    intensity, t, mol, peak_properties, isotopologue_lib
                )
                mol_i.append((mol, sum(intesity)))
                mol_peaks = list(zip(mz, intensity))
                scan_peaks.extend(mol_peaks)
                mol_scan_dict[mol]["ms1_scans"].append(i)
                # if len(mol_peaks) > 0:
                highest_peak = max(mol_peaks, key=lambda x: x[1])
                # else:
                #     highest_peak = (0, 0)
                mol_monoisotopic[mol] = {"mz": highest_peak[0], "i": highest_peak[1]}
        scan_peaks = sorted(scan_peaks, key=lambda x: x[1])
        if len(scan_peaks) > 0:
            mz, inten = zip(*scan_peaks)
        else:
            mz, inten = [], []
        s = Scan({"mz": np.array(mz), "i": np.array(inten), "id": i, "rt": t,})
        prec_scan_id = i
        i += 1
        scans.append((s, []))
        t += ms_rt_diff
        if t > gradient_length:
            break

        while len(mol_i) < 10:
            mol_i.append((None, 0))
        # for ms2 in range(10):
        for mol, _intensity in sorted(mol_i, key=lambda x: x[1], reverse=True)[:10]:
            mol_plus = f"{mol}"
            if mol is None:
                # add empty scan
                ms2_scan = Scan(
                    {
                        "mz": [],
                        "i": [],
                        "rt": t,
                        "id": i,
                        "precursor_mz": 0,
                        "precursor_i": 0,
                        "precursor_charge": 1,
                        "precursor_scan_id": prec_scan_id,
                    }
                )
                t += ms_rt_diff
                if t > gradient_length:
                    break
            elif (peak_properties[mol_plus]["scan_start_time"] <= t) and (
                (
                    peak_properties[mol_plus]["scan_start_time"]
                    + peak_properties[mol_plus]["peak_width"]
                )
                >= t
            ):
                mol_name = peak_properties[mol_plus]["trivial_name"]
                peaks = fragmentor.fragment(mol_name)
                frag_mz = peaks[:, 0]
                frag_i = peaks[:, 1]
                ms2_scan = Scan(
                    {
                        "mz": frag_mz,
                        "i": frag_i,
                        "rt": t,
                        "id": i,
                        "precursor_mz": mol_monoisotopic[mol]["mz"],
                        "precursor_i": mol_monoisotopic[mol]["i"],
                        "precursor_charge": 1,
                        "precursor_scan_id": prec_scan_id,
                    }
                )
                t += ms_rt_diff
                if t > gradient_length:
                    break
            # else:
            #     # print(t)
            #     ms2_scan = Scan(
            #         {
            #             "mz": [],
            #             "i": [],
            #             "rt": t,
            #             "id": i,
            #             "precursor_mz": None,
            #             "precursor_i": None,
            #             "precursor_charge": 1,
            #             "precursor_scan_id": prec_scan_id,
            #         }
            #     )
            #     t += ms_rt_diff
            # if t > gradient_length:
            #     break
            if mol is not None:
                mol_scan_dict[mol]["ms2_scans"].append(i)
            scans[-1][1].append(ms2_scan)
            i += 1
    return scans, mol_scan_dict


def generate_molecule_isotopologue_lib(
    peak_properties: Dict[str, dict],
    charges: List[int] = None,
    trivial_names: Dict[str, str] = None,
):
    """Summary.

    Args:
        molecules (TYPE): Description
    """
    if charges is None:
        charges = [1]
    if len(peak_properties) > 0:
        molecules = peak_properties.keys()
        lib = pyqms.IsotopologueLibrary(
            molecules=molecules,
            charges=charges,
            verbose=False,
            trivial_names=trivial_names,
        )
        reduced_lib = {}
        # for mol in lib:
        for mol in molecules:
            formula = lib.lookup["molecule to formula"][mol]
            data = lib[formula]["env"][(("N", "0.000"),)]
            # todo use correct charge states ==> hardcoded 1
            # try:
            #     triv_name = lib.lookup['molecule to trivial name'][mol]
            # except KeyError:
            #     triv_name = lib.lookup['molecule to trivial name'][f"+{mol}"]
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
    id_format_str = "controllerType=0 controllerNumber=1 scan={i}"
    with MzMLWriter(file) as writer:
        # Add default controlled vocabularies
        writer.controlled_vocabularies()
        # Open the run and spectrum list sections
        time_array = []
        intensity_array = []
        with writer.run(id="Simulated Run"):
            spectrum_count = len(scans) + sum([len(products) for _, products in scans])
            with writer.spectrum_list(count=spectrum_count):
                for scan, products in scans:
                    # Write Precursor scan
                    try:
                        index_of_max_i = np.argmax(scan.i)
                        max_i = scan.i[index_of_max_i]
                        mz_at_max_i = scan.mz[index_of_max_i]
                    except ValueError:
                        mz_at_max_i = 0
                        max_i = 0
                    spec_tic = sum(scan.i)
                    writer.write_spectrum(
                        scan.mz,
                        scan.i,
                        id=id_format_str.format(i=scan.id),
                        params=[
                            "MS1 Spectrum",
                            {"ms level": 1},
                            {"scan start time": scan.retention_time},
                            {"total ion current": spec_tic},
                            {"base peak m/z": mz_at_max_i},
                            {"base peak intensity": max_i},
                        ],
                    )
                    time_array.append(scan.retention_time)
                    intensity_array.append(spec_tic)
                    # Write MSn scans
                    for prod in products:
                        writer.write_spectrum(
                            prod.mz,
                            prod.i,
                            id=id_format_str.format(i=prod.id),
                            params=[
                                "MSn Spectrum",
                                {"ms level": 2},
                                {
                                    "scan start time": scan.retention_time,
                                    "unitName": "seconds",
                                },
                                {"total ion current": sum(prod.i)},
                            ],
                            # TOFIX adding precursor information makes psims crash
                            # Include precursor information
                            # precursor_information=prec_info
                            precursor_information={
                                "mz": prod.precursor_mz,
                                "intensity": prod.precursor_i,
                                "charge": prod.precursor_charge,
                                "scan_id": id_format_str.format(i=scan.id),
                                "spectrum_reference": id_format_str.format(i=scan.id),
                                "activation": ["HCD", {"collision energy": 25.0}],
                            },
                        )
            with writer.chromatogram_list(count=1):
                writer.write_chromatogram(
                    time_array,
                    intensity_array,
                    id="TIC",
                    chromatogram_type="total ion current",
                )
    return
