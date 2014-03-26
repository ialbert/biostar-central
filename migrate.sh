#!/bin/bash
#
# Performs an export from a biostar 1 site.
# (about 8 minutes for 100K posts)
#
#
# Usage: migrate filename
#


set -ue


FNAME=$1

REMOTE="biostar@biostars.org:~/data/$FNAME"
LOCAL=import/$FNAME

echo "*** fetching $REMOTE"

# Pull the latest file from remote, this is a postgres datadump
rsync -avz --rsh='ssh' $REMOTE $LOCAL

# Prepare the environment
source conf/default.env
source conf/migrate.env

# Drop the old database and import the new values
./biostar.sh pgdrop pgcreate pgimport $LOCAL

# Export the data into a directory
python -m main.bin.export -u -p -v -d /tmp/full
