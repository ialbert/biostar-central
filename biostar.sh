#!/bin/bash

# Environment variables must be set externally.
if [ -z "$BIOSTAR_HOME" ]; then
    echo "(!) environment variables not set."
    echo "(!) try: source conf/defaults.env"
    exit 1
fi

# Stop on errors or missing environment variables.
set -ue

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
    echo '  import   - imports a JSON data fixture'
    echo '  dump     - dumps the current database as a JSON data fixture'
    echo '  delete   - removes the sqlite database (sqlite specific)'
    echo '  run      - runs the development server'
    echo '  test     - runs all tests'
    echo '  env      - shows all customizable environment variables'
    echo ''
    echo "Use environment variables to customize settings. See the docs."
    echo ''
    echo "DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE"
    echo ''
fi

while (( "$#" )); do

    if [ "$1" = "run" ]; then
        echo "*** run the development server"
        $PYTHON $DJANGO_ADMIN runserver $BIOSTAR_HOSTNAME --settings=$DJANGO_SETTINGS_MODULE
    fi

    if [ "$1" = "delete" ]; then
        echo "*** deleting the sqlite database"
        $PYTHON $DJANGO_ADMIN delete_database --settings=$DJANGO_SETTINGS_MODULE
    fi

    if [ "$1" = "init" ]; then
        echo "*** initializing server on $BIOSTAR_HOSTNAME"

        $PYTHON $DJANGO_ADMIN test -v $VERBOSITY --settings=$DJANGO_SETTINGS_MODULE
        $PYTHON $DJANGO_ADMIN syncdb -v $VERBOSITY --noinput --settings=$DJANGO_SETTINGS_MODULE

        #$PYTHON_EXE $DJANGO_ADMIN migrate main.server --settings=$DJANGO_SETTINGS_MODULE
        #$PYTHON_EXE $DJANGO_ADMIN migrate djcelery --settings=$DJANGO_SETTINGS_MODULE
        #$PYTHON_EXE $DJANGO_ADMIN migrate kombu.transport.django --settings=$DJANGO_SETTINGS_MODULE
        #echo "*** collecting static files"
        $PYTHON $DJANGO_ADMIN collectstatic -v $VERBOSITY --noinput --settings=$DJANGO_SETTINGS_MODULE
    fi

    # Produce the environment variables recognized by Biostar.
    if [ "$1" = "test" ]; then
        echo "*** running all test"
        $PYTHON $DJANGO_ADMIN test -v $VERBOSITY --settings=$DJANGO_SETTINGS_MODULE
    fi

    # Produce the environment variables recognized by Biostar.
    if [ "$1" = "env" ]; then
        echo "*** Biostar specific environment variables"
        echo BIOSTAR_HOME=$BIOSTAR_HOME
        echo BIOSTAR_STATIC_ROOT=$BIOSTAR_STATIC_ROOT
        echo DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE
    fi

    if [ "$1" = "import_biostar1" ]; then
        echo "*** migrating from Biostar 1"
        echo "*** BIOSTAR_MIGRATE_DIR=$BIOSTAR_MIGRATE_DIR"
        $PYTHON $DJANGO_ADMIN import_biostar1 -u -p -x
    fi

shift
done