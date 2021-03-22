#!/bin/bash

# Default database backup  script.

cd /export/www/biostar-central/

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

export POSTGRES_HOST=/var/run/postgresql

# Activate the conda environemnt.
conda activate engine

# Stop on errors.
set -ue

USER=www

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

# Backup location
mkdir -p export/backup

# pg_dump the database
python manage.py tasks --action pg_dump --outdir export/backup  --user ${USER} --hourly

# bump a random post
# python manage.py tasks --action bump
