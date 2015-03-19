#!/bin/bash

# Environment variables must be set externally.
if [ -z "$BIOSTAR_HOME" ]; then
    echo "(!) Environment variables not set. See the README.md."
    echo "(!) Try: source run/sqlite.env"
    exit 1
fi

# Optionall override the python executable.
PYTHON=${PYTHON:=python}
VERBOSITY=${VERBOSITY:=1}

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
    echo '  init      - initializes the database'
    echo '  run       - runs the development server'
    echo "  index     - initializes the search index"
    echo '  test      - runs all tests'
    echo '  env       - shows all customizable environment variables'
    echo ' '
    #echo "  import    - imports the data fixture JSON_DATA_FIXTURE=$JSON_DATA_FIXTURE"
    #echo "  dump      - dumps data as JSON_DATA_FIXTURE=$JSON_DATA_FIXTURE"
    echo "  delete    - removes the sqlite database DATABASE_NAME=$DATABASE_NAME"
    echo ''
    echo "  pg_drop   - drops postgres DATABASE_NAME"
    echo "  pg_create - creates postgres DATABASE_NAME"
    echo "  pg_import file.gz - imports the gzipped file into postgres DATABASE_NAME"
    echo ''
    echo "Use environment variables to customize settings. See the docs."
    echo ' '
    echo "DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE"
    echo "DATABASE_NAME=$DATABASE_NAME"
    echo ''
fi

while (( "$#" )); do

     if [ "$1" = "init" ]; then
        echo "*** Initializing server on $BIOSTAR_HOSTNAME with $DJANGO_SETTINGS_MODULE"

        $PYTHON manage.py migrate --settings=$DJANGO_SETTINGS_MODULE

        #echo "*** Running all tests"
        #$PYTHON manage.py test --noinput -v $VERBOSITY --settings=$DJANGO_SETTINGS_MODULE
        #$PYTHON manage.py syncdb -v $VERBOSITY --noinput --settings=$DJANGO_SETTINGS_MODULE

        #$PYTHON manage.py migrate  biostar.apps.users --settings=$DJANGO_SETTINGS_MODULE
        #$PYTHON manage.py migrate  biostar.apps.posts --settings=$DJANGO_SETTINGS_MODULE
        #$PYTHON manage.py migrate  --settings=$DJANGO_SETTINGS_MODULE
        #$PYTHON manage.py initialize_site --settings=$DJANGO_SETTINGS_MODULE

        #$PYTHON manage.py collectstatic -v $VERBOSITY --noinput --settings=$DJANGO_SETTINGS_MODULE
    fi

	if [ "$1" = "quick" ]; then
		# Used during development only.
        echo "*** Quick init: delete import init loaddata"
		$PYTHON manage.py patch --delete_sqlite --settings=$DJANGO_SETTINGS_MODULE

 		# Load the database dump into sqlite.
 		sqlite3 $DATABASE_NAME < init/biostar2/biostar2-sqlite3.sql

		# Migrate the database.
        $PYTHON manage.py migrate --settings=$DJANGO_SETTINGS_MODULE

		# You will need to create this file - a dump of social authentication data - see docs.
		$PYTHON manage.py loaddata init/socialapp.json
    fi

    if [ "$1" = "run" ]; then
        echo "*** Run the development server with $DJANGO_SETTINGS_MODULE and DATABASE_NAME=$DATABASE_NAME"
        $PYTHON manage.py runserver $BIOSTAR_HOSTNAME --settings=$DJANGO_SETTINGS_MODULE
    fi

    if [ "$1" = "delete" ]; then
		# Deletes the sqlite database. Used during development.
        $PYTHON manage.py patch --delete_sqlite --settings=$DJANGO_SETTINGS_MODULE
    fi

 	# Produce the environment variables recognized by Biostar.
    if [ "$1" = "test" ]; then
        echo "*** Running all tests"
        $PYTHON manage.py test biostar3 --noinput --failfast -v $VERBOSITY --settings=run.sqlite
    fi

    if [ "$1" = "pg_drop" ]; then
        echo "*** Dropping the $DATABASE_NAME=$DATABASE_NAME!"
        dropdb -i $DATABASE_NAME
    fi

    if [ "$1" = "pg_create" ]; then
        # creates the PG database
        echo "*** Creating postgresql database DATABASE_NAME=$DATABASE_NAME"
        createdb $DATABASE_NAME -E utf8 --template template0
    fi

 	if [ "$1" = "import" ]; then
        echo "*** Importing data with DATA_IMPORT_COMMAND"
        eval $DATA_IMPORT_COMMAND
    fi

    if [ "$1" = "pg_import" ]; then
        echo "*** Importing into DATABASE_NAME=$DATABASE_NAME"
        gunzip -c $2 | psql $DATABASE_NAME
    fi

    if [ "$1" = "pg_dump" ]; then
        echo "*** Dumping the $DATABASE_NAME database."
        $PYTHON manage.py biostar_pg_dump -v $VERBOSITY --settings=$DJANGO_SETTINGS_MODULE
    fi

	if [ "$1" = "waitress" ]; then
        echo "*** Run a waitress server with $DJANGO_SETTINGS_MODULE and DATABASE_NAME=$DATABASE_NAME"
        waitress-serve --port=8080 biostar3.wsgi_whitenoise:application
    fi

   	if [ "$1" = "testdeploy" ]; then
        echo "*** deploys to the test site"
        fab -f conf/fabs/fabfile.py test_site pull restart
    fi

    # Produce the environment variables recognized by Biostar.
    if [ "$1" = "env" ]; then
        echo "*** Biostar specific environment variables"
        echo BIOSTAR_HOME=$BIOSTAR_HOME
        echo BIOSTAR_ADMIN_EMAIL=$BIOSTAR_ADMIN_EMAIL
        echo BIOSTAR_ADMIN_NAME=$BIOSTAR_ADMIN_NAME
        echo "-"
        echo DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE
        echo DATABASE_NAME=$DATABASE_NAME
        echo DEFAULT_FROM_EMAIL=$DEFAULT_FROM_EMAIL
    fi

    if [ "$1" = "dump" ]; then
        echo "*** Dumping json data into $JSON_DATA_FIXTURE"
        $PYTHON manage.py dumpdata users posts messages badges planet --settings=$DJANGO_SETTINGS_MODULE | gzip > $JSON_DATA_FIXTURE
    fi

    if [ "$1" = "index" ]; then
        echo "*** Indexing site content"
        $PYTHON manage.py rebuild_index --noinput --settings=$DJANGO_SETTINGS_MODULE
    fi

    if [ "$1" = "update_index" ]; then
        echo "*** Updating site index"
        $PYTHON manage.py update_index --age 1 --settings=$DJANGO_SETTINGS_MODULE
    fi

    if [ "$1" = "import_biostar1" ]; then
        echo "*** Migrating from Biostar 1"
        echo "*** BIOSTAR_MIGRATE_DIR=$BIOSTAR_MIGRATE_DIR"
        $PYTHON manage.py import_biostar1 -u -p -x
    fi


shift
done