#
# pyinfra backup.py file
#
from pyinfra.operations import files

from vars import *

files.directory(
    name="Create the directory",
    path=dirname(PG_DATA_DEST),
)

files.rsync(
    name="Copy over the latest backup file",
    src=PG_DATA_SRC,
    dest=PG_DATA_DEST,
    flags=['-a', ]
)
