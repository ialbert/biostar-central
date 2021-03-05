#!/bin/bash

# Default database backup  script.

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

export POSTGRES_HOST=/var/run/postgresql

# Activate the conda environemnt.
conda activate engine

USER=www

# Stop on errors.
set -ue

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

# Backup location
mkdir -p export/backup

# pg_dump the database
python manage.py tasks --action pg_dump --outdir export/backup  --user ${USER}


