import os
from pipeline.const import *

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
# To initialize test files run:
#
# make testfile
#

TEST_PROJECTS = [
    ("Sequencing run 1", "Lamar sequencing center"),
    ("Sequencing run 2", "Lamar sequencing center"),
    ("Sequencing run 3", "Lamar sequencing center"),
]

TEST_FASTQ_PATH = os.path.expandvars('tmp/1-SarriPal_S9_L001_R1_001.fastq.gz')
TEST_FILE_PATH = os.path.expandvars('tmp/sampleinfo.txt')
TEST_COLL_PATH = os.path.expandvars('tmp/data.tar.gz')
TEST_DATA = [
    ("Data Collection 1", "This file contains a collection of data", TEST_COLL_PATH, TAR_FASTQ_GZ),
    ("Data Collection 2", "This file contains a collection of data", TEST_COLL_PATH, TAR_FASTQ_GZ),
    ("Sample sheet 1", "This file contains a sample sheet describing the data in the directory", TEST_FILE_PATH,
     SAMPLE_TYPE),
    ("Sample sheet 2", "This file contains a sample sheet describing the data in the directory", TEST_FILE_PATH,
     SAMPLE_TYPE),
    ("Fastq file 1", "This is a fastq file", TEST_FASTQ_PATH, FASTQ_TYPE),
    ("Fastq file 2", "This is a fastq file", TEST_FASTQ_PATH, FASTQ_TYPE),

]

TEST_SPECS = [
    join(__THIS, '..', 'pipeline', 'qc', 'qc.hjson'),
    join(__THIS, '..', 'pipeline', 'fastqc', 'fastqc.hjson')
]


# User information

LOGIN_URL = "/login/"
