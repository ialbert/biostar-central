"""
Main test script is executed when running::

    biostar.sh test"

"""
import sys, os, django

if django.VERSION < (1, 3):
    print '*** Django version 1.3 or higher required.'
    print '*** Your version is %s' % '.'.join( map(str, django.VERSION))
    sys.exit()

from django.utils import unittest
from django.test.simple import DjangoTestSuiteRunner
from coverage import coverage
  
cov = coverage(include = [ 'main/*' ], omit=[ 'main/server/tests/*' ] )

COVERAGE = 0
if COVERAGE:
    cov.start()

# add our own testing suites
from main.server import html
from main.server.tests import test_models, test_site

def path(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))

# the directory that this file is located in
__CURR_DIR  = path(os.path.dirname(__file__))
REPORT_PATH = path(__CURR_DIR, '..', '..', '..', 'report')

class BiostarTest(DjangoTestSuiteRunner):

    def __init__(self, verbosity=1, interactive=True, failfast=True, **kwargs):
        super( BiostarTest, self ).__init__(verbosity, interactive, failfast, **kwargs)

    def run_tests(self, test_labels, extra_tests=None, **kwargs):
        # add new tests then delegate to supercalss
        
        extra_tests = [
            test_site.suite(),  
            html.suite(),
            test_models.suite(), 
        ]

        code = super( BiostarTest, self ).run_tests(test_labels, extra_tests, **kwargs)
        if COVERAGE and code == 0:
            cov.stop()
            # reporting if there are not errors
            cov.report()
            cov.html_report(directory=REPORT_PATH)
            #cov.xml_report()

