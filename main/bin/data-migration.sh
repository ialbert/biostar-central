#!/bin/bash

#
# migrates an existing datadump
# between schema changes using an intermediate checkout
#
#
# start this command in afresh checkout of biostar
#
# you will also need to have proper environments variables
#
# careful with the settings - this script drops databases!
#
set -ue

#
# make a copy of this file and customize as necessary
#
source conf/default

export PG_DBNAME=biostar-database
export PG_USER=biostar-user
export PG_PASSWD=biostar-password

PGDUMP=updated-biostar.sql

START_REV=0938eac
START_ENV=migconf/migrate_start.env

# this is currently ignored and the latest 
# revision is taken
END_REV=de9a834
END_ENV=migconf/migrate_end.env

#this is only necessary if the migration is repeated
#git branch -D start-migration

# check out the schema before the database migration
git checkout -b start-migration $START_REV

# initialize the environment
source $START_ENV

# initialize wit the old data then migrate
./biostar.sh pgdrop pgcreate pgimport

# get the most up to date version
git checkout master

# apply the new settings
source libs/$END_ENV

python manage.py syncdb
python manage.py migrate main.server 0001 --fake
python manage.py migrate main.server

echo "*** apply new ranking"
python -m main.bin.patch --reapply_ranks --update_domain

echo "*** dumping data to $PGDUMP"
./biostar.sh pgdump > $PGDUMP

echo "*** done"