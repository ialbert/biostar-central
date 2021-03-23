#!/bin/bash

# Stop on errors.
set -ue

# Postgres user.
USER=www

# Postgres database
DATABASE=biostardb

# Backup directory.
DIR=/export/www/biostar-central/export/backup

# Make the backup directory.
mkdir -p $DIR

# Build the name for the file
TSTAMP=`date +"daily-%Y-%m-%d"`

# Build the name of the backup file.
NAME=$DATABASE-$TSTAMP.gz

# Backup filename.
FILE=$DIR/$NAME

# Dump the database.
/usr/bin/pg_dump -x -O -b -U $USER  $DATABASE | gzip > $FILE

# Check the files the pattern finds.
#find $DIR/*.gz -mtime +4 -exec echo {} \;

# Dangerous command! Make sure it is correct!
#find $DIR/*.gz -mtime +4 -exec rm {} \;
