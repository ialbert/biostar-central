#!/bin/bash
set -ue

# verbosity level for commands 0=minimal, 2=maximal
VERBOSITY=1

# set a few default environment variables
BIOSTAR_SRC=`dirname $0`

# source directory to be added to import path
BIOSTAR_HOME=${BIOSTAR_HOME:-"$BIOSTAR_SRC/main"}

# set the hostname
BIOSTAR_HOSTNAME=${BIOSTAR_HOSTNAME:-"0.0.0.0:8080"}

# django settings module
export DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE:-"settings"}

# add to the python path
PYTHONPATH=${PYTHONPATH:-""}

# the migration path and limit
export MIGRATE_PATH=${MIGRATE_PATH:-"import/se2"}
export MIGRATE_LIMIT=${MIGRATE_LIMIT:-"100"}

# the fixture to dump/load data from
export FIXTURE=${FIXTURE:-"import/datadump.json"}
export FIXTURE_GZ=$FIXTURE.gz

# the DJANGO_SETTINGS_MODULE needs to be in the python import path
export PYTHONPATH=$PYTHONPATH:$BIOSTAR_HOME   

# add the library files to the pythonpath
export PYTHONPATH=$PYTHONPATH:libs/:libs/libraries.zip

# setting up the python
export PYTHON_EXE=python
export DJANGO_ADMIN=main/manage.py

echo ""
echo "Settings:"
echo "*** BIOSTAR_HOME=$BIOSTAR_HOME"
echo "*** BIOSTAR_HOSTNAME=$BIOSTAR_HOSTNAME"
echo "*** DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE"
echo "*** PYTHONPATH=$PYTHONPATH"

if [ $# == 0 ]; then
	echo ''
	echo Usage:
	echo '  $ run.sh <command>'
	echo ''
	echo 'Multiple commands may be used in the same line:'
	echo '  $ run.sh init import run'
	echo ''
	echo 'Commands:'
	echo '  init     - initializes the database'
	echo '  populate - populates the system with test data'
	echo '  delete   - removes everything from BioStar'
	echo '  run      - runs server'
	echo '  import   - imports from a StackExchange data export'
	echo '  test     - runs all tests'

fi

while (( "$#" )); do

	if [ "$1" = "delete" ]; then
                # hardcore reset that is needed if tables need to reinitialized
		echo "*** deleting sqlite"
		rm -f $BIOSTAR_HOME/db/biostar.db
	fi

        if [ "$1" = "drop" ]; then
		echo '*** dropping postgresql'
		dropdb biostar-test-database
                echo '*** creating postgresql'
                createdb biostar-test-database
	fi
        
	if [ "$1" = "flush" ]; then
                echo "*** flushing the database"
		$PYTHON_EXE $DJANGO_ADMIN flush --noinput --settings=$DJANGO_SETTINGS_MODULE
		
	fi

	if [ "$1" = "init" ]; then
		echo "*** initializing server on $BIOSTAR_HOSTNAME"
		$PYTHON_EXE $DJANGO_ADMIN syncdb -v $VERBOSITY --noinput --settings=$DJANGO_SETTINGS_MODULE

                echo "*** collecting static files"
                $PYTHON_EXE $DJANGO_ADMIN collectstatic -v $VERBOSITY --noinput --settings=$DJANGO_SETTINGS_MODULE
	fi

	if [ "$1" = "import" ]; then
		echo "*** imports data from $FIXTURE_GZ"
		$PYTHON_EXE $DJANGO_ADMIN loaddata $FIXTURE_GZ --settings=$DJANGO_SETTINGS_MODULE
		echo "*** indexing post content"
		$PYTHON_EXE -m main.server.search  --settings=$DJANGO_SETTINGS_MODULE
	fi

	if [ "$1" = "run" ]; then
		echo "*** running the webserver on $BIOSTAR_HOSTNAME"
		$PYTHON_EXE $DJANGO_ADMIN runserver $BIOSTAR_HOSTNAME --settings=$DJANGO_SETTINGS_MODULE
	fi

	if [ "$1" = "test" ]; then
		echo "*** running the tests"
		#$PYTHON_EXE $DJANGO_ADMIN test server --settings=$DJANGO_SETTINGS_MODULE --failfast
		$PYTHON_EXE $DJANGO_ADMIN test server --settings=$DJANGO_SETTINGS_MODULE --failfast
	fi

	if [ "$1" = "dump" ]; then		
		echo "*** dumping data to $FIXTURE"
		$PYTHON_EXE $DJANGO_ADMIN dumpdata auth.User server --settings=$DJANGO_SETTINGS_MODULE | gzip > $FIXTURE_GZ
	fi

	if [ "$1" = "migrate" ]; then
		echo "*** migrating data to a new datadump"
		source conf/memory.sh
		echo "*** DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE"
		echo "*** MIGRATE_PATH=$MIGRATE_PATH, $MIGRATE_LIMIT"
		echo "*** FIXTURE=$FIXTURE_GZ"
		$PYTHON_EXE -m main.migrate -o $FIXTURE --path $MIGRATE_PATH --limit $MIGRATE_LIMIT
		gzip -f $FIXTURE
		echo "*** dumped data to $FIXTURE_GZ"
	fi

	if [ "$1" = "index" ]; then		
		echo "*** indexing all post content"
		$PYTHON_EXE -m main.server.search
	fi

shift
done

