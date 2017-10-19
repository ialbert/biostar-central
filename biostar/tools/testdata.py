import os

from . import const


def join(*args):
    abspath = os.path.abspath(os.path.join(*args))
    return abspath

# Data needs to be stored relative to this directory
__THIS = os.path.dirname(__file__)
__ROOT = join(__THIS, "..", "..")

# Used only during testing.
TEST_PROJECTS = [
    ("Sequencing run 3", "Lamar sequencing center"),
    ("Sequencing run 2", "Lamar sequencing center"),
    ("Sequencing run 1", "Lamar sequencing center"),
]

TEST_FASTQ = join(__ROOT, 'tmp/1-SarriPal_S9_L001_R1_001.fastq.gz')
TEST_SAMPLE = join(__ROOT, 'tmp/sampleinfo.txt')
TEST_COLLECTION = join(__ROOT, 'tmp/data.tar.gz')

if not os.path.isfile(TEST_FASTQ):
    print("*** Missing test data. Run: make testdata")

TEST_DATA = [
    ("Data Collection 1.tar.gz", "A test collection of data", TEST_COLLECTION, const.TAR_FASTQ_GZ),
    ("Data Collection 2.tar.gz", "A test collection of data", TEST_COLLECTION, const.TAR_FASTQ_GZ),
    ("Sample sheet 1.txt", "A test sample sheet.", TEST_SAMPLE, const.SAMPLE_TYPE),
    ("Sample sheet 2.txt", "A test sample sheet", TEST_SAMPLE, const.SAMPLE_TYPE),
    ("FASTQ file 2.fq", "This is a fastq file", TEST_FASTQ, const.FASTQ_TYPE),
    ("FASTQ file 2.fq", "This is a fastq file", TEST_FASTQ, const.FASTQ_TYPE),
]
