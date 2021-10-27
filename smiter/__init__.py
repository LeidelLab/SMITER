"""Top-level package for SMITER."""

__author__ = """Manuel Kösters"""
__email__ = "manuel.koesters@dcb.unibe.ch"
__version__ = "0.1.0"

import os
import sys
import tempfile

from loguru import logger

import smiter.synthetic_mzml

config = {
    "handlers": [
        {"sink": sys.stdout, "level": "INFO"},
        {
            "sink": os.path.join(tempfile.gettempdir(), "smiter.log"),
            "compression": "gz",
            "rotation": "100 MB",
            "level": "DEBUG",
            "backtrace": True,
        },
    ],
}
logger.configure(**config)
