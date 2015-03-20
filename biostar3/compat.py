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

try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError

def strip(text):
    return text.strip()