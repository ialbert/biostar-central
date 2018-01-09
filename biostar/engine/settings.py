import hjson
import os
from biostar.settings import *

DATA_TYPES_FILE = hjson.load(open(os.path.abspath("./conf/datatypes.hjson")))

DATA_TYPES = [
              (name, DATA_TYPES_FILE[name].get("symbol", ''), DATA_TYPES_FILE[name].get("help",""))
              for name in DATA_TYPES_FILE
             ]



