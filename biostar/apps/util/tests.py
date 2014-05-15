"""
Utility related tests.

These will execute when you run "manage.py test".
"""
from __future__ import print_function, unicode_literals, absolute_import, division

import bleach, logging
from django.conf import settings
from biostar.apps.util import html

from django.test import TestCase

logging.disable(logging.INFO)
# The pattern that matches the user link.

class UtilTest(TestCase):

    def test_bleach(self):
        "Testing html cleaning"
        eq = self.assertEqual

        inp = '''http://www.psu.edu'''
        exp = '''<a href="http://www.psu.edu" rel="nofollow">http://www.psu.edu</a>'''
        got = bleach.linkify(inp)
        #eq(got, exp)

        inp = '''http://%s/u/123''' % settings.SITE_DOMAIN
        exp = '''<a href="http://www.psu.edu" rel="nofollow">http://www.psu.edu</a>'''
        got = bleach.linkify(inp)
        #eq(got, exp)

        #print (bleach.DEFAULT_CALLBACKS)




