#!/bin/bash

BATCH_SIZE=5000

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

# Stop on errors.
set -ue

export POSTGRES_HOST=/var/run/postgresql

python manage.py index --clear_spam