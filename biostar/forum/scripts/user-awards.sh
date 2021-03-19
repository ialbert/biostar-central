#!/bin/bash


cd /export/www/biostar-central/

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

export POSTGRES_HOST=/var/run/postgresql

# Activate the conda environemnt.
conda activate engine

# Stop on errors.
set -ue

LIMIT=100

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

python manage.py tasks --action award --limit ${LIMIT}