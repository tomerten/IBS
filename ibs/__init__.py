# -*- coding: utf-8 -*-

"""
Package ibs
=======================================

Top-level package for ibs.
"""

import ibs.cli_runode

__version__ = "0.0.0"


import IBSLib as ibslib

from .cli_runode import _MODEL_MAP, plot


def print_model_map() -> None:
    """Print the IBS model map."""
    print(_MODEL_MAP)
