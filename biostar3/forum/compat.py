"""
Compatibility with Python 3.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from io import open

try:
    import itertools.imap as map
    import itertools.ifilter as filter
    import itertools.izip as zip
except ImportError:
    pass


def strip(text):
    return text.strip()