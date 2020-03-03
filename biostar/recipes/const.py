import os

__THIS = os.path.dirname(__file__)


def join(*args):
    return os.path.abspath(os.path.join(*args))


KNOWN_TEXT_EXTENSIONS = (
    ".fasta", ".fq", ".fastq", ".sam", ".bed", ".gff", ".fa",
    ".r", ".gtf", ".gb", ".tsv", ".csv", ".rep", ".sh", ".md", ".txt"
)


# Stop words used in search
STOPWORDS = (
    'a,able,about,across,after,all,almost,also,am,among,an,'
    'and,any,are,as,at,be,because,been,but,by,can,cannot,'
    'could,dear,did,do,does,either,else,ever,every,for,from,'
    'get,got,had,has,have,he,her,hers,him,his,how,however,i,'
    'if,in,into,is,it,its,just,least,let,like,likely,may,me,'
    'might,most,must,my,neither,no,nor,not,of,off,often,on,'
    'only,or,other,our,own,rather,said,say,says,she,should,'
    'since,so,some,than,that,the,their,them,then,there,'
    'these,they,this,tis,to,too,twas,us,wants,was,we,were,'
    'what,when,where,which,while,who,whom,why,will,with,'
    'would,yet,you,your').split(',')

REDIRECT_FIELD_NAME = 'next'

COPIED_DATA = "data"
COPIED_RESULTS = "results"
COPIED_RECIPES = "recipes"
COPIED_FILES = "files"
CLONED_RECIPES = "clone"

CLIPBOARD_CONTENTS = [COPIED_DATA, COPIED_RESULTS, COPIED_RECIPES, COPIED_FILES,
                      CLONED_RECIPES]

# Paste as data objects
DATA_PASTE_TARGETS = ','.join([COPIED_DATA, COPIED_RESULTS, COPIED_FILES])

# Paste as recipe objects
RECIPE_PASTE_TARGETS = ','.join([COPIED_RECIPES, CLONED_RECIPES])

KNOWN_TEXT_EXTENSIONS = set(KNOWN_TEXT_EXTENSIONS)

# The maximum length in characters for a typical name and text field.
MAX_NAME_LEN = 256
MAX_FIELD_LEN = 1024
MAX_TEXT_LEN = 10000
MAX_LOG_LEN = 20 * MAX_TEXT_LEN

# Map a file extension to a biostar-engine datatype.
EXT_TO_TYPE = dict(

    fa='FASTA',
    fasta='FASTA',
    txt='TXT',
    fq='FQ',
    gtf='GTF',
    gb='GB',
    tsv='TSV',
    csv='CSV',
    rep='REP',
    md='MD',
    sam='SAM',
    bed='BED',
    gff='GFF',
    r='R',
    fastq='FASTQ'

)

RADIO = "RADIO"
DROPDOWN = "DROPDOWN"
INTEGER = "INTEGER"
UPLOAD = "UPLOAD"
FLOAT = "FLOAT"
CHECKBOX = "CHECKBOX"
TEXTBOX = "TEXTBOX"
SQL = "SQL"
