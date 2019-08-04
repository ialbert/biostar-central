#!/bin/bash

# Database migration script.

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

# Activate the conda environemnt.
conda activate engine

# Stop on errors.
set -ue

export DJANGO_SETTINGS_MODULE=themes.bioconductor.settings

# Migrate the server.
python manage.py migrate

# Collect static files
python manage.py collectstatic --noinput
