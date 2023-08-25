import os

def path(p):
    p = os.path.abspath(os.path.expanduser(p))
    return p

# Source code directory
SOURCE_DIR = path('~/app/biostar-central')

PG_DATA_NAME = 'biostardb-daily-2023-08-08.gz'

# The name of the remote file
PG_DATA_SRC = f'www@test.biostars.org:/home/www/data/biostars/{PG_DATA_NAME}'

# The directory for the backups
PG_DATA_DIR = path(f'{SOURCE_DIR}/export/backup')

# The name of the backup
PG_DATA_DEST = path(f'{PG_DATA_DIR}/{PG_DATA_NAME}')

# The name of the local file
PG_DATABASE_NAME = "biostar"
PG_ROLE = 'www'
PG_ROLE_PASSWORD = 'test'

PSQL_USER = 'ialbert'
PSQL_PASSWORD = 'test'
PSQL_HOST = 'localhost'



def dirname(path):
    return os.path.dirname(path)
