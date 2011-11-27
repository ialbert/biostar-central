#!/bin/bash
set -ue

# verbosity level for commands 0=minimal, 2=maximal
VERBOSITY=1

# set a few default environment variables
BIOSTAR_SRC=`dirname $0`

# source directory to be added to import path
BIOSTAR_HOME=${BIOSTAR_HOME:-"$BIOSTAR_SRC/main"}

# this will contain the collected static files
BIOSTAR_EXPORT=${BIOSTAR_EXPORT:-"$BIOSTAR_SRC/export"}

# set the hostname
BIOSTAR_HOSTNAME=${BIOSTAR_HOSTNAME:-"0.0.0.0:8080"}

# django settings module
DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE:-"settings"}

# add to the python path
PYTHONPATH=${PYTHONPATH:-""}

# the fixture to dump/load data from
export FIXTURE=import/datadump.json.gz

# the DJANGO_SETTINGS_MODULE needs to be in the python import path
export PYTHONPATH=$PYTHONPATH:$BIOSTAR_HOME   

# add the library files to the pythonpath
export PYTHONPATH=$PYTHONPATH:libs/

# setting up the python
export PYTHON_EXE=python
export DJANGO_ADMIN=main/manage.py

echo ""
echo "Settings:"
echo "*** BIOSTAR_HOME=$BIOSTAR_HOME"
echo "*** BIOSTAR_EXPORT=$BIOSTAR_EXPORT"
echo "*** BIOSTAR_HOSTNAME=$BIOSTAR_HOSTNAME"
echo "*** DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE"

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
		echo '*** deleting all data'
		# list all commands that need to be run
		cmds[0]="rm -f $BIOSTAR_HOME/db/biostar.db"
		for cmd in "${cmds[@]}"
			do
				echo "*** executing: $cmd"
				`$cmd`
			done
	fi

	if [ "$1" = "init" ]; then
		echo "*** initializing server on $BIOSTAR_HOSTNAME"
		$PYTHON_EXE $DJANGO_ADMIN syncdb -v $VERBOSITY --noinput --settings=$DJANGO_SETTINGS_MODULE
                 
                echo "*** collecting static files to $BIOSTAR_EXPORT"
                mkdir -p $BIOSTAR_EXPORT
                $PYTHON_EXE $DJANGO_ADMIN collectstatic -v $VERBOSITY --noinput --settings=$DJANGO_SETTINGS_MODULE
	fi

	if [ "$1" = "populate" ]; then
		echo "*** populating server with: $FIXTURE"
		$PYTHON_EXE $DJANGO_ADMIN loaddata $FIXTURE --settings=$DJANGO_SETTINGS_MODULE
	fi

	if [ "$1" = "run" ]; then
		echo "*** running the webserver on $BIOSTAR_HOSTNAME"
		$PYTHON_EXE $DJANGO_ADMIN runserver $BIOSTAR_HOSTNAME --settings=$DJANGO_SETTINGS_MODULE
	fi

	if [ "$1" = "test" ]; then
		echo "*** running the tests"
		$PYTHON_EXE $DJANGO_ADMIN test server --settings=$DJANGO_SETTINGS_MODULE --failfast
	fi

	if [ "$1" = "flush" ]; then
		echo "*** flushing data"
		$PYTHON_EXE $DJANGO_ADMIN flush --settings=$DJANGO_SETTINGS_MODULE > $FIXTURE
	fi

	if [ "$1" = "dump" ]; then		
		echo "*** dumping data to $FIXTURE"
		$PYTHON_EXE $DJANGO_ADMIN dumpdata auth.User server --settings=$DJANGO_SETTINGS_MODULE | gzip > $FIXTURE
	fi

	if [ "$1" = "import" ]; then
		echo "*** importing the data into the main database"
		#rm -f $BIOSTAR_HOME/db/biostar.db
		#cp $BIOSTAR_HOME/db/test.db $BIOSTAR_HOME/db/biostar.db
		$PYTHON_EXE import/migrate.py --path import/se0 --limit 100
		#$PYTHON_EXE import/migrate.py --path import/se2 --limit 500
		#$PYTHON_EXE import/migrate.py --path import/se2
	fi

shift
done

