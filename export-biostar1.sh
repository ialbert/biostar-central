#!/bin/bash
#
# Performs an export from a biostar 1 site.
# (about 8 minutes for 100K posts)
#
#
# Usage: migrate filename
#

set -ue

FNAME=latest.sql.gz

REMOTE_SQL="biostar@biostars.org:~/data/$FNAME"
LOCAL_SQL=import/$FNAME
MIGRATE_DIR=/tmp/migrate

echo "*** fetching $REMOTE_SQL to $LOCAL_SQL"

# Pull the latest file from remote, this is a postgres datadump
rsync -avz --rsh='ssh' $REMOTE_SQL $LOCAL_SQL

# Prepare the environment
source conf/default.env
source conf/migrate.env

# Drop the old database and import the new values
./biostar.sh pgdrop pgcreate pgimport $LOCAL_SQL

# Export the data into a directory
python -m main.bin.export -u -p -v -d $MIGRATE_DIR -n 100
