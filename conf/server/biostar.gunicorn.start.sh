#!/bin/bash
set -ue

# This will set the environment variables.
source /home/www/sites/biostar-central/live/deploy.env

PID="/home/www/sites/biostar-central/live/biostar.gunicorn.pid"

# The user and group the unicorn process will run as.
NUM_WORKERS=3

# Where to bind.
BIND="unix:/tmp/biostar.sock"
#BIND="localhost:8080"

# The WSGI module that starts the process.
DJANGO_WSGI_MODULE='biostar.wsgi'

# The gunicorn instance to run.
GUNICORN="gunicorn"

# How many requests to serve.
MAX_REQUESTS=1000

# The name of the application.
NAME="biostar_app"

# --log-level=debug \

exec $GUNICORN ${DJANGO_WSGI_MODULE}:application \
  --name $NAME \
  --workers $NUM_WORKERS \
  --max-requests $MAX_REQUESTS\
  --bind $BIND
  --pid $PID\
