"""
Main test script is executed when running::

    biostar.sh test"

"""
import sys, os, unittest, django

if django.VERSION < (1, 3):
    print '*** Django version 1.3 or higher required.'
    print '*** Your version is %s' % '.'.join( map(str, django.VERSION))
    sys.exit()

from django.test.simple import DjangoTestSuiteRunner

try:
    from coverage import coverage
except ImportError:
    coverage = None

# add our own testing suites
from biostar.server import html
from biostar.tests import functional, access

class BiostarTest(DjangoTestSuiteRunner):

    def __init__(self, verbosity=1, interactive=True, failfast=True, **kwargs):
        super( BiostarTest, self ).__init__(verbosity, interactive, failfast, **kwargs)

    def run_tests(self, test_labels, extra_tests=None, **kwargs):
        # add new tests then delegate to supercalss
        extra_tests = [  
            access.suite(), functional.suite(), html.suite(),
        ]

        if coverage:
            cov = coverage(include = ['biostar/*'], omit=['biostar/libs/*', 'biostar/tests/*'] )
            cov.start()
            super( BiostarTest, self ).run_tests(test_labels, extra_tests, **kwargs)
            cov.stop()
            cov.report()
            cov.html_report()
            cov.xml_report()
        else:
            super( BiostarTest, self ).run_tests(test_labels, extra_tests, **kwargs)


