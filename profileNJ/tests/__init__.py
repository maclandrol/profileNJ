# This file is part of profileNJ testing module

__author__ = "Emmanuel Noutahi"

import unittest
import sys
import os

dirname = os.path.dirname(__file__)

test_suite = unittest.defaultTestLoader.discover(__name__)
