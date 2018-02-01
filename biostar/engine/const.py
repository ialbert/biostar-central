import os
from biostar.tools.const import *

__THIS = os.path.dirname(__file__)


def join(*args):
    return os.path.abspath(os.path.join(*args))



