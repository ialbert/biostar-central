"""
Compatibility with Python 3.
"""

from io import open

try:
    import itertools.imap as map
    import itertools.ifilter as filter
    import itertools.izip as zip
except ImportError:
    pass
