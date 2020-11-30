#!/bin/bash

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

export POSTGRES_HOST=/var/run/postgresql

# Activate the conda environment.
conda activate engine

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

# Check for user awards every hour
python manage.py awards --awards_users