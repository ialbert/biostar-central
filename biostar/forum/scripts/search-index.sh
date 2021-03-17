#!/bin/bash

BATCH_SIZE=5000

cd /export/www/biostar-central/

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

export POSTGRES_HOST=/var/run/postgresql

# Activate the conda environemnt.
conda activate engine

# Stop on errors.
set -ue

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

# Set the log
# export ENGINE_LOG_LEVEL=DEBUG

# Add BATCH_SIZE posts to search index
python manage.py index --size ${BATCH_SIZE}