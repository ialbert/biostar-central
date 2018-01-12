import hjson
import os
from biostar.settings import *


def expand(dict):
    return (dict.get("name", ''), dict.get("symbol", ""), dict.get("help", ""))

DATA_TYPES_FILE = os.path.abspath("./conf/datatypes.hjson")

TYPE_DICT = hjson.load(open(DATA_TYPES_FILE))

DATA_TYPES = [ expand(n) for n in TYPE_DICT ]





