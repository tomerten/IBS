# -*- coding: utf-8 -*-

"""Tests for ibs package."""

import ibs


def test_hello_noargs():
    """Test for ibs.hello()."""
    s = ibs.hello()
    assert s=="Hello world"


def test_hello_me():
    """Test for ibs.hello('me')."""
    s = ibs.hello('me')
    assert s=="Hello me"
    
    
# ==============================================================================
# The code below is for debugging a particular test in eclipse/pydev.
# (otherwise all tests are normally run with pytest)
# Make sure that you run this code with the project directory as CWD, and
# that the source directory is on the path
# ==============================================================================
if __name__ == "__main__":
    the_test_you_want_to_debug = test_hello_noargs

    print("__main__ running", the_test_you_want_to_debug)
    the_test_you_want_to_debug()
    print('-*# finished #*-')
    
# eof