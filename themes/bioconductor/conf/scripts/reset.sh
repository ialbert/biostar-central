#!/bin/bash

# Default database reset script.

# Stop on errors.
set -ue

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

# Activate the correct enviroment.
conda activate engine

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=themes.bioconductor.settings

#python manage.py flush --noinput
make reset
