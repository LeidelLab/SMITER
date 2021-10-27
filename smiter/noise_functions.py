"""Callables for injection noise into scans.

Upon calling the callabe, a list/np.array of mz and intensities should be returned.
Arguments should be passed via *args and **kwargs
"""
import math
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple

import numpy as np
import pyqms
from loguru import logger

import smiter
from smiter.lib import calc_mz

# from smiter.synthetic_mzml import Scan
# import smiter.synthetic_mzml


class AbstractNoiseInjector(ABC):
    """Summary."""

    def __init__(self, *args, **kwargs):
        """Initialize noise injector."""
        pass  # pragma: no cover

    @abstractmethod
    def inject_noise(self, scan, *args, **kwargs):
        """Main noise injection method.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description

        """
        pass  # pragma: no cover

    @abstractmethod
    def _ms1_noise(self, scan, *args, **kwargs):
        pass  # pragma: no cover

    @abstractmethod
    def _msn_noise(self, scan, *args, **kwargs):
        pass  # pragma: no cover


class GaussNoiseInjector(AbstractNoiseInjector):
    def __init__(self, *args, **kwargs):
        # np.random.seed(1312)
        logger.info("Initialize GaussNoiseInjector")
        self.args = args
        self.kwargs = kwargs

    def inject_noise(self, scan, *args, **kwargs):
        """Main noise injection method.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description

        """
        self.kwargs.update(kwargs)
        # self.args += args
        if scan.ms_level == 1:
            scan = self._ms1_noise(scan, *self.args, **self.kwargs)
        elif scan.ms_level > 1:
            scan = self._msn_noise(scan, *args, **kwargs)
        return scan

    def _ms1_noise(self, scan, *args, **kwargs):
        """Generate ms1 noise.

        Args:
            scan (Scan): Description
            *args: Description
            **kwargs: Description
        """
        mz_noise = self._generate_mz_noise(scan, *args, **kwargs)
        intensity_noise = self._generate_intensity_noise(scan, *args, **kwargs)
        scan.mz += mz_noise
        scan.i += intensity_noise
        scan.i[scan.i < 0] = 0
        return scan

    def _msn_noise(self, scan, *args, **kwargs):
        """Generate msn noise.

        Args:
            scan (Scan): Description
            *args: Description
            **kwargs: Description
        """
        mz_noise = self._generate_mz_noise(scan, *args, **kwargs)
        intensity_noise = self._generate_intensity_noise(scan, *args, **kwargs)
        scan.mz += mz_noise
        scan.i += intensity_noise
        scan.i[scan.i < 0] = 0
        # remove some mz values
        dropout = kwargs.get("dropout", 0.1)
        dropout_mask = [
            False if c < dropout else True for c in np.random.random(len(scan.mz))
        ]
        scan.mz = scan.mz[dropout_mask]
        scan.i = scan.i[dropout_mask]
        return scan

    def _generate_mz_noise(self, scan, *args, **kwargs):
        """Generate noise for mz_array.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description
        """
        ppm_offset = kwargs.get("ppm_offset", 0)
        ppm_var = kwargs.get("ppm_var", 1)
        noise_level = np.random.normal(0, ppm_var * 1e-6, len(scan.mz))
        # noise_level = np.zeros(len(scan.mz))
        noise = scan.mz * noise_level
        # noise = np.zeros(len(scan.mz))
        return noise

    def _generate_intensity_noise(self, scan, *args, **kwargs):
        """Generate intensity noise.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description
        """
        noise = np.random.normal(
            # np.array(np.zeros(len(scan.i))),
            np.array(scan.i),
            np.array(scan.i) * kwargs.get("variance", 0.02),
            len(scan.i),
        )
        return noise - scan.i


class UniformNoiseInjector(AbstractNoiseInjector):
    def __init__(self, *args, **kwargs):
        # np.random.seed(1312)
        logger.info("Initialize UniformNoiseInjector")
        self.args = args
        self.kwargs = kwargs

    def inject_noise(self, scan, *args, **kwargs):
        """Main noise injection method.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description

        """
        self.kwargs.update(kwargs)
        # self.args += args
        if scan.ms_level == 1:
            scan = self._ms1_noise(scan, *self.args, **self.kwargs)
        elif scan.ms_level > 1:
            scan = self._msn_noise(scan, *self.args, **self.kwargs)
        return scan

    def _ms1_noise(self, scan, *args, **kwargs):
        """Generate ms1 noise.

        Args:
            scan (Scan): Description
            *args: Description
            **kwargs: Description
        """
        mz_noise = self._generate_mz_noise(scan, *args, **kwargs)
        intensity_noise = self._generate_intensity_noise(scan, *args, **kwargs)
        scan.mz += mz_noise
        scan.i += intensity_noise
        scan.i[scan.i < 0] = 0
        return scan

    def _msn_noise(self, scan, *args, **kwargs):
        """Generate msn noise.

        Args:
            scan (Scan): Description
            *args: Description
            **kwargs: Description
        """
        mz_noise = self._generate_mz_noise(scan, *args, **kwargs)
        intensity_noise = self._generate_intensity_noise(scan, *args, **kwargs)
        scan.mz += mz_noise
        scan.i += intensity_noise
        scan.i[scan.i < 0] = 0
        # scan.mz[scan.mz < 0] = 0
        dropout = kwargs.get("dropout", 0.1)
        dropout_mask = [
            False if c < dropout else True for c in np.random.random(len(scan.mz))
        ]
        scan.mz = scan.mz[dropout_mask]
        scan.i = scan.i[dropout_mask]
        return scan

    def _generate_mz_noise(self, scan, *args, **kwargs):
        """Generate noise for mz_array.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description
        """
        # noise_level = np.random.normal(ppm_offset * 5e-6, ppm_var * 5e-6, len(scan.mz))
        # get scaling from kwargs
        noise = np.random.uniform(
            (scan.mz * kwargs.get("ppm_noise", 5e-6)) * -1,
            scan.mz * kwargs.get("ppm_noise", 5e-6),
        )
        # noise = scan.mz * noise_level
        return noise

    def _generate_intensity_noise(self, scan, *args, **kwargs):
        """Generate intensity noise.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description
        """
        # get scaling from kwargs
        noise = np.random.uniform(
            (scan.i * kwargs.get("intensity_noise", 0.2)) * -1,
            scan.i * kwargs.get("intensity_noise", 0.2),
        )
        return noise


class PPMShiftInjector(AbstractNoiseInjector):
    def __init__(self, *args, **kwargs):
        # np.random.seed(1312)
        logger.info("Initialize PPMShiftInjector")
        self.args = args
        self.kwargs = kwargs

    def inject_noise(self, scan, *args, **kwargs):
        """Main noise injection method.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description

        """
        self.kwargs.update(kwargs)
        # self.args += args
        if scan.ms_level == 1:
            scan = self._ms1_noise(scan, *self.args, **self.kwargs)
        elif scan.ms_level > 1:
            scan = self._msn_noise(scan, *self.args, **self.kwargs)
        return scan

    def _ms1_noise(self, scan, *args, **kwargs):
        """Generate ms1 noise.

        Args:
            scan (Scan): Description
            *args: Description
            **kwargs: Description
        """
        mz_noise = self._generate_mz_noise(scan, *args, **kwargs)
        intensity_noise = self._generate_intensity_noise(scan, *args, **kwargs)
        scan.mz += mz_noise
        scan.i += intensity_noise
        scan.i[scan.i < 0] = 0
        return scan

    def _msn_noise(self, scan, *args, **kwargs):
        """Generate msn noise.

        Args:
            scan (Scan): Description
            *args: Description
            **kwargs: Description
        """
        mz_noise = self._generate_mz_noise(scan, *args, **kwargs)
        intensity_noise = self._generate_intensity_noise(scan, *args, **kwargs)
        scan.mz += mz_noise
        scan.i += intensity_noise
        scan.i[scan.i < 0] = 0
        # scan.mz[scan.mz < 0] = 0
        dropout = kwargs.get("dropout", 0.1)
        dropout_mask = [
            False if c < dropout else True for c in np.random.random(len(scan.mz))
        ]
        scan.mz = scan.mz[dropout_mask]
        scan.i = scan.i[dropout_mask]
        return scan

    def _generate_mz_noise(self, scan, *args, **kwargs):
        """Generate noise for mz_array.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description
        """
        # noise_level = np.random.normal(ppm_offset * 5e-6, ppm_var * 5e-6, len(scan.mz))
        # get scaling from kwargs
        offset = self.kwargs["offset"]
        sigma = self.kwargs["sigma"]
        noise = np.random.normal(offset, sigma, len(scan.mz))
        noise = noise * scan.mz
        return noise

    def _generate_intensity_noise(self, scan, *args, **kwargs):
        """Generate intensity noise.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description
        """
        # get scaling from kwargs
        noise = [0 for x in range(len(scan.i))]
        return noise


# TODO add white noise params
# TODO should be called MSpireNoiseInjector
class JamssNoiseInjector(AbstractNoiseInjector):
    def __init__(self, *args, **kwargs):
        """Noise injector based on this paper:
            https://academic.oup.com/bioinformatics/article/31/5/791/318378

        Args:
            **kwargs: fine tuning for mz and itensity sigma
                - m: scale mz noise sigma by this value
                - y: controll falling of mz sigma with higher intensities with this
                - a: scale intensity sigma by this
                - b: controll falling of intensity sigma with higher intensities with this
                - c: constant to add to intensity sigma
        """
        # np.random.seed(1312)
        logger.info("Initialize JamssNoiseInjector")
        self.args = args
        self.kwargs = kwargs

    def inject_noise(self, scan, *args, **kwargs):
        self.kwargs.update(kwargs)
        if scan.ms_level == 1:
            scan = self._ms1_noise(scan, *self.args, **self.kwargs)
        elif scan.ms_level > 1:
            scan = self._msn_noise(scan, *self.args, **self.kwargs)
        return scan

    def _add_white_noise(self, scan):
        n = int(np.random.uniform(100, 500))
        if sum(scan.i) == 0:
            total_tic = 5e6
        else:
            total_tic = sum(scan.i)
            total_tic = 5e6
        # logger.info(f"Total TIC: {sum(scan.i)}")
        white_noise_i = np.random.uniform(1, 100, n)
        max_perc_noise = 0.75
        min_perc_noise = 0.5
        white_noise_i = (
            white_noise_i
            / white_noise_i.sum()
            * np.random.uniform(min_perc_noise, max_perc_noise)
            * total_tic
        )
        # logger.info(f'Add white noise intensity: {white_noise_i}')
        # scale = np.random.uniform(1e5, 1e8) / n
        # white_noise_i *= scale
        # logger.info(f"Noise: {sum(white_noise_i)}")
        white_noise_mz = np.random.uniform(0, 1200, n)
        new_i = np.concatenate((scan.i, white_noise_i))
        new_mz = np.concatenate((scan.mz, white_noise_mz))
        # logger.info(f'Add new white noise {new_i}')
        sort = new_mz.argsort()
        scan.mz = new_mz[sort]
        scan.i = new_i[sort]
        return scan

    def _ms1_noise(self, scan, *args, **kwargs):
        # i_b = sum(scan.i)
        scan = self._add_white_noise(scan)
        mz_noise = self._generate_mz_noise(scan, *args, **kwargs)
        # logger.debug(f"MS1 mz noise: {mz_noise}")
        intensity_noise = self._generate_intensity_noise(scan, *args, **kwargs)
        scan.mz += mz_noise
        scan.i += intensity_noise

        # i_a = sum(scan.i)
        # print(f'S2I: {i_b/i_a}')
        scan.i[scan.i < 0] = 0
        return scan

    def _msn_noise(self, scan, *args, **kwargs):
        mz_noise = self._generate_mz_noise(scan, *args, **kwargs)
        # logger.debug(f"MSn mz noise: {mz_noise}")
        intensity_noise = self._generate_intensity_noise(scan, *args, **kwargs)
        scan.mz += mz_noise
        scan.i += intensity_noise
        scan.i[scan.i < 0] = 0
        # scan.mz[scan.mz < 0] = 0
        dropout = kwargs.get("dropout", 0.1)
        # logger.debug(f"Dropout: {dropout}")
        dropout_mask = [
            False if c < dropout else True for c in np.random.random(len(scan.mz))
        ]
        scan.mz = scan.mz[dropout_mask]
        scan.i = scan.i[dropout_mask]
        scan = self._add_white_noise(scan)
        return scan

    def _generate_mz_noise(self, scan, *args, **kwargs):
        try:
            norm_int = scan.i / max(scan.i) * 100
        except ValueError:
            norm_int = np.array([])
        # TODO make variable parameters
        # TODO Check mspire paper again for variables
        # TODO pass data structure with max_i in elution profile and mzs calculate noise to intensity/max_i_in_profile
        m = 0.001701
        # m = 0.1701
        y = 0.543
        y = 0.2
        sigma = m * norm_int ** (-y)
        noise = np.random.normal(loc=scan.mz, scale=sigma)
        noise -= scan.mz
        noise_mock = np.zeros(len(scan.mz))
        return noise

    def _generate_intensity_noise(self, scan, *args, **kwargs):
        try:
            norm_int = scan.i / max(scan.i) * 100
        except ValueError:
            norm_int = np.array([])
        # TODO make variable parameters
        # TODO Check mspire paper again for variables
        # TODO pass data structure with max_i in elution profile and mzs calculate noise to intensity/max_i_in_profile
        m = 0.001701
        m = 10.34
        c = 0.00712
        d = 0.12
        sigma = m * (1 - math.e ** (-c * norm_int)) + d
        noise = np.random.normal(loc=0, scale=sigma)
        if len(scan.i) > 0:
            noise *= max(scan.i) / 100
        else:
            noise = np.array([])
        # TODO add white noise points with noise = x * tic (0 < x < 1)
        return noise
