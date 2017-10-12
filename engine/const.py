import os

__THIS = os.path.dirname(__file__)


def join(*args):
    return os.path.abspath(os.path.join(*args))


HOME_ICON = "home"
PROJECT_LIST_ICON = "database"
PROJECT_ICON = "archive"
DATA_LIST_ICON = "file text"
DATA_ICON = "file"
ANALYSIS_LIST_ICON = "settings"
ANALYSIS_ICON = "setting"
RESULT_LIST_ICON = "bar chart"
RESULT_ICON = "line chart"
FIELD_VISIBLE = "visible"
FIELD_ORIGIN = "origin"
PROJECT_ORIGIN = "PROJECT"

#
# Run to initialize test files.
#
# make testfile
#
TEST_FILE = os.path.expandvars('tmp/sampleinfo.txt')
TEST_COLL = os.path.expandvars('tmp/data.tar.gz')
TEST_DATA = [
    ("Compressed data directory", "This file contains a collection of data", TEST_COLL),
    ("Sample sheet", "This file contains a sample sheet describing the data in the directory", TEST_FILE),
]

TEST_SPECS = [
    join(__THIS, '..', 'pipeline', 'qc', 'qc_spec.hjson')
]
