"""Callables for injection noise into scans.

Upon calling the callabe, a list/np.array of mz and intensities should be returned.
Arguments should be passed via *args and **kwargs
"""
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple

import numpy as np
import pyqms

import smiter
from smiter.lib import calc_mz

# from smiter.synthetic_mzml import Scan
# import smiter.synthetic_mzml

easter_egg = """

                   .-^-.
                 .'=^=^='.
                /=^=^=^=^=\

        .-~-.  :^= HAPPY =^;
      .'~~*~~'.|^ EASTER! ^|
     /~~*~~~*~~\^=^=^=^=^=^:
    :~*~~~*~~~*~;\.-*))`*-,/
    |~~~*~~~*~~|/*  ((*   *'.
    :~*~~~*~~~*|   *))  *   *\

     \~~*~~~*~~| *  ((*   *  /
      `.~~*~~.' \  *))  *  .'
        `~~~`    '-.((*_.-'
"""


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
        print(easter_egg)
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
        self.args += args
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
        # breakpoint()
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
        print(easter_egg)
        self.args = args
        self.kwargs = kwargs

    def inject_noise(self, scan, *args, **kwargs):
        """Main noise injection method.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description

        """
        # breakpoint()
        self.kwargs.update(kwargs)
        self.args += args
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
        print(intensity_noise)
        scan.mz += mz_noise
        # breakpoint()
        scan.i += intensity_noise

        scan.i[scan.i < 0] = 0
        # scan.mz[scan.mz < 0] = 0
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
            (scan.mz * kwargs.get('ppm_noise', 5e-6)) * -1,
            scan.mz * kwargs.get('ppm_noise', 5e-6),
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
            (scan.i * kwargs.get('intensity_noise', 0.2)) * -1,
            scan.i * kwargs.get('intensity_noise', 0.2),
        )
        return noise
