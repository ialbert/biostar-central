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

PGDUMP=updated-biostar.sql

START_REV=0938eac
START_ENV=migconf/migrate_start.env

END_REV=de9a834
END_ENV=migconf/migrate_end.env

# check out the schema before the database migration
git checkout -b start-migration $START_REV

# initialize the environment
source $START_ENV

# initialize wit the old data then migrate
./biostar.sh pgdrop pgcreate pgimport

git checkout -b end-migration $END_REV

source $END_ENV

python manage.py syncdb
python manage.py migrate main.server --fake

echo "*** apply new ranking"
python -m main.bin.patch --reapply_ranks --update_domain

echo "*** dumping data to $PGDUMP"
./biostar.sh pgdump > $PGDUMP

echo "*** done"