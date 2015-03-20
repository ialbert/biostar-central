"""
Compatibility with Python 3.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from io import open

try:
    # Try the Python 2 versions first.
    import itertools.imap as map
    import itertools.ifilter as filter
    import itertools.izip as zip
except ImportError:
    # These all work correctly in Python 3.
    pass

try:
    # Try the Python 3 versions first.
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
except ImportError:
    # Use the Python 2 versions.
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError


def strip(text):
    return text.strip()

from io import StringIO
from django.utils.encoding import smart_text