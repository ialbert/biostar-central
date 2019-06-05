#!/bin/bash

# Default database reset script.

# Activate the correct enviroment.
source /home/www/miniconda3/envs/engine/bin/activate engine

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

#python manage.py flush --noinput
make reset
