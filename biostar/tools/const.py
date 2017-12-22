# Valid interface types
RADIO = "RADIO"
DROPDOWN = "DROPDOWN"
INTEGER = "INTEGER"
UPLOAD = "UPLOAD"
FLOAT = "FLOAT"
CHECKBOX = "CHECKBOX"
TEXTBOX = "TEXTBOX"

# Data types
GENERIC_TYPE = 1
SAMPLE_TYPE = 2
COLLECTION_TYPE = 3
TAR_FASTQ_GZ = 4
FASTQ_TYPE = 5
FASTA_TYPE = 6
REFERENCE_TYPE = 7
LAMAR_SAMPLE_SHEET = 8
ACCESSION_MAP = 9
ACCESSION_LIST = 10
GTF = 11

DATA_TUPLES = [
    # Numeric, Symbol, Readable English
    (GENERIC_TYPE, 'GENERIC', "Generic"),
    (SAMPLE_TYPE, 'SAMPLE_TYPE', "Sample Sheet"),
    (COLLECTION_TYPE, 'COLLECTION', "Collection"),
    (TAR_FASTQ_GZ, 'TAR_FASTQ_GZ', "TAR FASTQ GZ"),
    (FASTQ_TYPE, 'FASTQ', "FASTQ"),
    (FASTA_TYPE, 'FASTA', "FASTA"),
    (REFERENCE_TYPE, 'REFERENCE', "Reference"),
    (LAMAR_SAMPLE_SHEET, 'LAMAR_SAMPLE_SHEET', "Lamar Sample Sheet"),
    (ACCESSION_MAP, 'ACCESSION_MAP', "Accession Map"),
    (ACCESSION_LIST, 'ACCESSION_LIST', "Accession List"),
    (GTF, 'GTF', "GTF Annotations"),

]

# Data type in plain English.
DATA_TYPES = dict(
    (a, c) for (a, b, c) in DATA_TUPLES
)

# Data type to symbol.
DATA_TYPE_SYMBOLS = dict(
    (b, a) for (a, b, c) in DATA_TUPLES
)

# Sequence subtype.
SEQUENCE_TYPES = [FASTA_TYPE, FASTQ_TYPE]
COLLECTION_TYPES = [TAR_FASTQ_GZ]
