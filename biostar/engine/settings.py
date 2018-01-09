import hjson
from biostar.settings import *

DATA_TYPES_FILE = "conf/datatypes.hjson"

DATA_TYPES = [
              (name, name.get("symbol", ""), name.get("help",""))
              for name in hjson.load(DATA_TYPES_FILE)
             ]



