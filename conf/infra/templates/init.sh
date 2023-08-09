#/bin/bash

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

# Activate the conda environment.
conda activate engine

# Stop on errors.
set -uex

# Set the configuration module.
export DJANGO_SETTINGS_MODULE={{DJANGO_SETTINGS_MODULE}}

# Migrate the server.
python manage.py migrate

# Collect static files
python manage.py collectstatic --noinput
