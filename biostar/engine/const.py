import os

__THIS = os.path.dirname(__file__)

def join(*args):
    return os.path.abspath(os.path.join(*args))

KNOWN_TEXT_EXTENSIONS = (
    ".fasta", ".fq", ".fastq", ".sam", ".bed", ".gff", ".fa",
    ".r", ".gtf", ".gb", ".tsv", ".csv", ".rep", ".sh", ".md", ".txt"
)

REDIRECT_FIELD_NAME = 'next'

DATA_CLIPBOARD = "data_clipboard"
RESULTS_CLIPBOARD = "results_clipboard"
RECIPE_CLIPBOARD = "recipe_clipboard"

KNOWN_TEXT_EXTENSIONS = set(KNOWN_TEXT_EXTENSIONS)

# Map a file extension to a biostar-engine datatype.
EXT_TO_TYPE = dict(

            fa = 'FASTA',
            fasta = 'FASTA',
            txt = 'TXT',
            fq = 'FQ',
            gtf = 'GTF',
            gb = 'GB',
            tsv = 'TSV',
            csv = 'CSV',
            rep = 'REP',
            md = 'MD',
            sam = 'SAM',
            bed = 'BED',
            gff = 'GFF',
            r = 'R',
            fastq = 'FASTQ'

                )



RADIO = "RADIO"
DROPDOWN = "DROPDOWN"
INTEGER = "INTEGER"
UPLOAD = "UPLOAD"
FLOAT = "FLOAT"
CHECKBOX = "CHECKBOX"
TEXTBOX = "TEXTBOX"
