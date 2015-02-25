#!/bin/bash
set -ue

# This is required so that the default configuration file works.
source /home/www/sites/biostar-central/live/deploy.env

# Location of the log file
LOGFILE=/home/www/sites/biostar-central/live/logs/%n-%i.log

# The concurrency level
NUM_WORKERS=2

# How many tasks per child process
MAX_TASK=1000

# The name of the application.
APP="biostar"

# The logging level
LOGLEVEL=info

# The gunicorn instance to run.
CELERY="/home/www/.virtualenvs/biostar/bin/celery"

echo "starting celery worker with DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE"

exec $CELERY -A $APP worker -l info --maxtasksperchild $MAX_TASK --concurrency $NUM_WORKERS -f $LOGFILE

