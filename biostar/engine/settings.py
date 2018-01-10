import hjson
import os
from biostar.settings import *

DATA_TYPES_FILE = os.path.abspath("./conf/datatypes.hjson")

TYPE_DICT = hjson.load(open(DATA_TYPES_FILE))

DATA_TYPES = [
              (n, TYPE_DICT[n].get("symbol", ''), TYPE_DICT[n].get("help","")) for n in TYPE_DICT
             ]




