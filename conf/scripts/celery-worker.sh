#!/bin/bash
set -ue

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

export POSTGRES_HOST=/var/run/postgresql

# Activate the conda environemnt.
conda activate engine

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

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

# Celery instance to run
CELERY=" ~/miniconda3/envs/engine/bin/celery"

echo "starting celery worker with DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE"

exec $CELERY -A $APP worker -l info --maxtasksperchild $MAX_TASK --concurrency $NUM_WORKERS -f $LOGFILE

