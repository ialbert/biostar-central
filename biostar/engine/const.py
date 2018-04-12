import os

__THIS = os.path.dirname(__file__)

def join(*args):
    return os.path.abspath(os.path.join(*args))

KNOWN_EXTENSIONS = (
    ".fasta", ".fq", ".fastq", ".sam", ".bed", ".gff", ".fa",
    ".r", ".gtf", ".gb", ".tsv", ".csv", ".rep", ".sh", ".md", ".txt"
)

KNOWN_EXTENSIONS = set(KNOWN_EXTENSIONS)


REMAPPING = dict(
                fa = 'fasta'
            )




RADIO = "RADIO"
DROPDOWN = "DROPDOWN"
INTEGER = "INTEGER"
UPLOAD = "UPLOAD"
FLOAT = "FLOAT"
CHECKBOX = "CHECKBOX"
TEXTBOX = "TEXTBOX"
