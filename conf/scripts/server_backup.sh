#!/bin/bash

# Default database backup  script.

# Activate the correct enviroment.
source /home/www/miniconda3/envs/engine/bin/activate engine

export DJANGO_SETTINGS_MODULE=conf.site.site_settings

DUMP_FILE=export/database/recipe-backup.json
BACKUP_DUMP_FILE=export/database/recipe_backup_`date +'%s'`.json

python manage.py dumpdata --exclude contenttypes > $DUMP_FILE

cp -f $DUMP_FILE $BACKUP_DUMP_FILE

# Produce a datadump count as a reminder.
ls -1 export/database/*.json | wc -l
~
