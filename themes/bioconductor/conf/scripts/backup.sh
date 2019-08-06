#!/bin/bash

# Default database backup  script.

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

# Activate the conda environemnt.
conda activate engine

# Stop on errors.
set -ue

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=themes.bioconductor.settings
export DATABASE_NAME=bioconductor.db

# Backup location
mkdir -p export/backup

# Generate backup file location.
TSTAMP=`date +'%Y-%m-%d-%H-%m'`
BACKUP="export/backup/data-${TSTAMP}.json"

# Dump the data in the desired format.
python manage.py dumpdata --exclude contenttypes > $BACKUP

# Show the backup filename.
echo "$BACKUP"

