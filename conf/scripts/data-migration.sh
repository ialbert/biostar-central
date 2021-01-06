#!/bin/bash

# Database migration script that  migrates from old biostar to biostar next.

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

# Activate the conda environemnt.
conda activate engine

# Stop on errors.
set -ue

# Comment out on local machine
export POSTGRES_HOST=/var/run/postgresql

# Set the location for the postgres sql dump
OLD_DUMP_DIR=export/sql/

mkdir -p ${OLD_DUMP_DIR}

# Get the database dump for old database (postgres sql dump).
#rsync -avz www@test.biostars.org:/home/www/biostar-engine/export/sql/biostar-database-2.3.0-hourly-00.sql.gz ${OLD_DUMP_DIR}

#OLD_DATABASE_DUMP=${1:=OLD_DUMP_DIR/biostar-database-2.3.0-hourly-00.sql.gz}
OLD_DATABASE_DUMP=${OLD_DUMP_DIR}/biostar-database-2.3.0-hourly-00.sql.gz

# How many items to load
LIMIT=1000

# Set the old database.
export OLD_DATABASE=old_biostar_db

# Set the new database
export NEW_DATABASE=biostar-database

# Set the configuration module.
TRANSFER_SETTINGS_MODULE=biostar.transfer.settings

#------------------------------------------------------------------------------------------

# Drop the old biostar database if exists.
dropdb --if-exists ${OLD_DATABASE}
# Drop the new biostar next database if exists.
dropdb --if-exists ${NEW_DATABASE}

# Create a database for old data dump
createdb ${OLD_DATABASE}
# Create the new biostar next database.
createdb ${NEW_DATABASE}

#  Migrate newly created biostar next database
python manage.py migrate --settings ${TRANSFER_SETTINGS_MODULE}

# Dump data into old biostar database
cat ${OLD_DATABASE_DUMP} | gunzip -c | psql -d ${OLD_DATABASE}

# Transfer the old database into new database.
python manage.py transfer --limit $LIMIT --settings ${TRANSFER_SETTINGS_MODULE}

