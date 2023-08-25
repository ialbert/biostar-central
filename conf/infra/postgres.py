from pyinfra.operations import postgresql, files

from vars import *

# Initialize the database.
postgresql.database(
    database=PG_DATABASE_NAME,
    present=True,
    owner=None,
    encoding=None,
    lc_collate=None,
    lc_ctype=None,
    tablespace=None,
    template='',
    connection_limit=None,
    psql_user=PSQL_USER,
    psql_password=PSQL_PASSWORD,
    psql_host=PSQL_HOST,
)

# Create the directory for the backup.
files.directory(
    name=f"Create the directory: {PG_DATA_DIR}",
    path=PG_DATA_DIR,
)


# Copy over the latest backup file.
files.rsync(
    name="Copy over the latest backup file",
    src=PG_DATA_SRC,
    dest=PG_DATA_DEST,
    flags=['-a', ]
)

'''
# Import the data dump.
postgresql.load(
    name="Import the data dump",
    src=PG_DATA_DEST,
    database=PG_DATABASE_NAME,
    psql_user=PSQL_USER, psql_password=PSQL_PASSWORD,
)
'''
