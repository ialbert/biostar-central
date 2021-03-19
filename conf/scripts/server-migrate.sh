#!/bin/bash

# Database migration script.

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

cd /export/www/biostar-central/

# Activate the conda environemnt.
conda activate engine

# Stop on errors.
set -ue

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

export POSTGRES_HOST=/var/run/postgresql

# Migrate the server.
python manage.py migrate

# Collect static files
python manage.py collectstatic --noinput
