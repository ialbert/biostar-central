#!/bin/bash

# Stop on errors or missing environment variables.
set -ue

# Setting up the default environment variables.
# Override them in the calling environment.

# Hostname for the development server.
BIOSTAR_HOSTNAME=${BIOSTAR_HOSTNAME:="localhost"}

# The settings files to be used.
DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE:="biostar.settings.devel"}

# The python executable to invoke.
PYTHON=python

# The level of verbosity for django commands.
VERBOSITY=${VERBOSITY:="1"}

# The django manager to run.
DJANGO_ADMIN=manage.py

if [ $# == 0 ]; then
    echo ''
    echo 'Usage:'
    echo ''
    echo "  $ $(basename $0) <command>"
    echo ''
    echo 'Multiple commands may be used on the same line:'
    echo ''
    echo "  $ $(basename $0) init import run"
    echo ''
    echo 'Commands:'
    echo ''
    echo '  init     - initializes the database'
    echo '  import   - imports a data fixture'
    echo '  dump     - dumps the current database as a data fixture'
    echo '  delete   - removes the sqlite database (sqlite specific)'
    echo '  run      - runs server'
    echo '  test     - runs all tests'
    echo '  env      - shows all customizable environment variables'
    echo ''
    echo "Use environment variables to customize settings. See docs"
    echo ''
    echo "DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE"
    echo ''
fi


while (( "$#" )); do

    if [ "$1" = "delete" ]; then
    	echo "*** deleting the sqlite database"
        #$PYTHON $DJANGO_ADMIN delete_database --settings=$DJANGO_SETTINGS_MODULE
    fi

    if [ "$1" = "init" ]; then
        echo "*** initializing server on $BIOSTAR_HOSTNAME"
        $PYTHON $DJANGO_ADMIN syncdb -v $VERBOSITY --noinput --settings=$DJANGO_SETTINGS_MODULE
        #$PYTHON_EXE $DJANGO_ADMIN migrate main.server --settings=$DJANGO_SETTINGS_MODULE
        #$PYTHON_EXE $DJANGO_ADMIN migrate djcelery --settings=$DJANGO_SETTINGS_MODULE
        #$PYTHON_EXE $DJANGO_ADMIN migrate kombu.transport.django --settings=$DJANGO_SETTINGS_MODULE
        #echo "*** collecting static files"
        $PYTHON $DJANGO_ADMIN collectstatic -v $VERBOSITY --noinput --settings=$DJANGO_SETTINGS_MODULE
    fi

shift
done