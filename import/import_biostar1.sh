#!/bin/bash
set -ue

LOCAL_SQL=import/biostar-migration.sql.gz
REMOTE_SQL='www@test.biostars.org:~/data'
IMPORT_DIR=/tmp/migrate
PG_USER=ialbert

source live/deploy.env

# Initialize the database
./biostar.sh pg_drop pg_create init

# Import the biostar1 datadump
python manage.py import_biostar1 -u -p -x -d $IMPORT_DIR

# Reset the sql sequences
python manage.py sqlfix --reset

# If this was the main site then the migration is done.

# If the migration took place on another site then one needs to dump the data
# and transfer it to the new site
/usr/bin/pg_dump -Fp -x -O -b -U $PG_USER $DATABASE_NAME | gzip > $LOCAL_SQL

scp $LOCAL_SQL $REMOTE_SQL