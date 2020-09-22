#!/bin/bash

# Default database reset script.

# Activate the correct enviroment.
source /home/www/miniconda3/envs/engine/bin/activate engine

export POSTGRES_HOST=/var/run/postgresql

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

# Set the site domain.
SITE_DOMAIN=${1}
export SITE_DOMAIN=${SITE_DOMAIN}

#python manage.py flush --noinput
make reset
