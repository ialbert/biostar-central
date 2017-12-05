# Valid interface types
RADIO = "RADIO"
DROPDOWN = "DROPDOWN"
INTEGER = "INTEGER"
UPLOAD = "UPLOAD"
FLOAT = "FLOAT"
CHECKBOX = "CHECKBOX"
TEXTBOX = "TEXTBOX"

# Data types
GENERIC_TYPE = 0
SAMPLE_TYPE = 1
COLLECTION_TYPE = 2
TAR_FASTQ_GZ = 3
FASTQ_TYPE = 4
FASTA_TYPE = 5
REFERENCE_TYPE = 6
LAMAR_SAMPLE_SHEET = 7
ACCESSION_MAP = 8
ACCESSION_LIST = 9


DATA_TUPLES = [
    # Numeric, Symbol, Readable English
    (GENERIC_TYPE, 'GENERIC', "Generic Type"),
    (SAMPLE_TYPE, 'SAMPLE_TYPE', "Sample Sheet"),
    (COLLECTION_TYPE, 'COLLECTION', "Collection"),
    (TAR_FASTQ_GZ, 'TAR_FASTQ_GZ', "TAR FASTQ GZ"),
    (FASTQ_TYPE, 'FASTQ', "FASTQ Type"),
    (FASTA_TYPE, 'FASTA', "FASTA Type"),
    (REFERENCE_TYPE, 'REFERENCE', "Reference"),
    (LAMAR_SAMPLE_SHEET, 'LAMAR_SAMPLE_SHEET', "Lamar Sample Sheet"),
    (ACCESSION_MAP, 'ACCESSION_MAP', "Accession Map"),
    (ACCESSION_LIST, 'ACCESSION_LIST', "Accession List"),

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
