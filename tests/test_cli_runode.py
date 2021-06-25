#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for `ibs.cli_runode ` CLI."""

import os

from click.testing import CliRunner
from ibs.cli_runode import main

from tests.test_cpp_ODE import THIS_DIR

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def test_main_all():
    runner = CliRunner()
    result = runner.invoke(main, [os.path.join(THIS_DIR, "sim_input_test.json")])
    assert result.exit_code == 0


def test_main_single():
    runner = CliRunner()
    result = runner.invoke(main, [os.path.join(THIS_DIR, "sim_input_test_single.json")])
    assert result.exit_code == 0


# ==============================================================================
# The code below is for debugging a particular test in eclipse/pydev.
# (normally all tests are run with pytest)
# ==============================================================================
if __name__ == "__main__":
    the_test_you_want_to_debug = test_main_all

    print(f"__main__ running {the_test_you_want_to_debug}")
    the_test_you_want_to_debug()
    print("-*# finished #*-")
# eof
