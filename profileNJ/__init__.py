# This file is part of profileNJ.

__author__ = "Emmanuel Noutahi"

__project__ = 'ProfileNJ'
__version__ = '20161018'

VERSION = __project__ + '-' + __version__


def test(verbosity=1):
    import unittest
    from .tests import test_suite

    unittest.TextTestRunner(verbosity=verbosity).run(test_suite)
