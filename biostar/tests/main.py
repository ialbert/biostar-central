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

# add our own testing suites
from biostar.tests import functional, access

COMPUTE_COVERAGE = False
if COMPUTE_COVERAGE:
    from coverage import coverage
    cov = coverage()
    cov.start()

class BiostarTest(DjangoTestSuiteRunner):

    def __init__(self, verbosity=1, interactive=True, failfast=True, **kwargs):
        super( BiostarTest, self ).__init__(verbosity, interactive, failfast, **kwargs)

    def run_tests(self, test_labels, extra_tests=None, **kwargs):
        # add new tests then delegate to supercalss
        extra_tests = [  access.suite(), functional.suite() ]
        super( BiostarTest, self ).run_tests(test_labels, extra_tests, **kwargs)




