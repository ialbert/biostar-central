#!/bin/bash
set -ue

# This is required so that the default configuration file works.
source /home/www/sites/biostar-central/live/deploy.env

# Setting the various access logs.
ACCESS_LOG=/home/www/sites/biostar-central/live/logs/celery-access.log
ERROR_LOG=/home/www/sites/biostar-central/live/logs/celery-error.log

# The user and group the unicorn process will run as.
NUM_WORKERS=3

# The name of the application.
NAME="biostar"

# The gunicorn instance to run.
CELERY="/home/www/.virtualenvs/biostar/bin/celery"

echo "celery starting with DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE"

exec $CELERY celery -A biostar.server.tasks worker -l info

