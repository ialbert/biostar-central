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

NEW_PGDUMP=updated-biostar.sql

OLD_REV=0938eac
OLD_ENV=mig/migrate-old.env

NEW_REV=de9a834
NEW_ENV=mig/migrate.env

# initialize the environment
source $OLD_ENV

# roll back to the last valid format
git checkout $OLD_REV

# initialize wit the old data then migrate
./biostar.sh pgdrop pgcreate pgimport

git checkout $NEW_REV
source $NEW_ENV

python manage.py syncdb
python manage.py migrate main.server --fake

echo "*** apply new ranking"
python -m main.bin.patch --reapply_ranks --update_domain

echo "*** dumping data to $NEW_PGDUMP"
./biostar.sh pgdump > $NEW_PGDUMP
echo "*** done"