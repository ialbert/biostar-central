#!/bin/bash

# This is an example database migrations script.

source activate engine

python manage.py migrate --settings conf.site.site_settings

