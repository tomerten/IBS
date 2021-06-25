# -*- coding: utf-8 -*-

"""
Package ibs
=======================================

Top-level package for ibs.
"""

import ibs.cli_runode
__version__ = "0.0.0"


def hello(who="world"):
    """'Hello world' method.

    :param str who: whom to say hello to
    :returns: a string
    """
    result = "Hello " + who
    return result


import IBSLib as ibslib
