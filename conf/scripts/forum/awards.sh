#!/bin/bash

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

export POSTGRES_HOST=/var/run/postgresql

# Activate the conda environemnt.
conda activate engine

LIMIT=100
# Stop on errors.
set -ue

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

# Add 5000 posts to search index every 3 minutes
python manage.py tasks --action award --limit ${LIMIT}