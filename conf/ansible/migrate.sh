#!/bin/bash

# This is an example database migrations script.

# Activate the correct enviroment.
source /home/www/miniconda3/envs/engine/bin/activate engine

# Set the configuration module.
DJANGO_SETTINGS_MODULE=conf.site.site_settings

#python manage.py flush --noinput
python manage.py migrate
python manage.py collectstatic --noinput -v 0
python manage.py project --json initial/tutorial/tutorial-project.hjson --privacy public --jobs
python manage.py project --root ../biostar-recipes --json projects/cookbook/cookbook-project.hjson --privacy public --jobs


# Add users as needed
# python manage.py add_user initial/initial-users.csv
# python manage.py add_access initial/initial-access.csv
