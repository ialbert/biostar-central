#!/bin/bash
set -ue

# set a few default environment variables
BIOSTAR_SRC=`dirname $0`

# source directory to be added to import path
BIOSTAR_HOME=${BIOSTAR_HOME:-"$BIOSTAR_SRC/main"}

# set the hostname
BIOSTAR_HOSTNAME=${BIOSTAR_HOSTNAME:-"0.0.0.0:8080"}

# django settings module
DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE:-"settings"}

# add to the python path
PYTHONPATH=${PYTHONPATH:-""}

# the fixture to dump/load data from
export FIXTURE=import/datadump.json

# the DJANGO_SETTINGS_MODULE needs to be in the python import path
export PYTHONPATH=$PYTHONPATH:$BIOSTAR_HOME   

# setting up the python
export PYTHON_EXE=python
export DJANGO_ADMIN=main/manage.py

echo "*** BIOSTAR_HOME=$BIOSTAR_HOME"
echo "*** BIOSTAR_HOSTNAME=$BIOSTAR_HOSTNAME"
echo "*** DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE"

if [ $# == 0 ]; then
	echo ''
	echo Usage:
	echo '  $ biostar.sh <command>'
	echo ''
	echo 'Multiple commands may be used in the same line:'
	echo '  $ biostar.sh init test run'
	echo ''
	echo 'Commands:'
	echo '  init - initializes the database'
	echo '  test - runs all tests'
	echo '  populate -populates the system with test data'
	echo '  delete - removes everything from BioStar'
	echo '  run - runs server'
fi

while (( "$#" )); do

	if [ "$1" = "delete" ]; then
		echo '*** Deleting all data'
		# list all commands that need to be run
		cmds[0]="rm -f $BIOSTAR_HOME/db/biostar.db"

		for cmd in "${cmds[@]}"
			do
				echo $cmd
				`$cmd`
			done
	fi

	if [ "$1" = "init" ]; then
		echo "*** Initializing server on $BIOSTAR_HOSTNAME"
		$PYTHON_EXE $DJANGO_ADMIN syncdb --noinput --settings=$DJANGO_SETTINGS_MODULE
	fi

	if [ "$1" = "populate" ]; then
		echo "*** Populating server with: $FIXTURE"
		$PYTHON_EXE $DJANGO_ADMIN loaddata $FIXTURE --settings=$DJANGO_SETTINGS_MODULE
	fi

	if [ "$1" = "run" ]; then
		echo "*** Running the webserver on $BIOSTAR_HOSTNAME"
		$PYTHON_EXE $DJANGO_ADMIN runserver $BIOSTAR_HOSTNAME --settings=$DJANGO_SETTINGS_MODULE
	fi

	if [ "$1" = "test" ]; then
		echo "*** Running the tests"
		$PYTHON_EXE $DJANGO_ADMIN test server --settings=$DJANGO_SETTINGS_MODULE --failfast
	fi

	if [ "$1" = "flush" ]; then
		echo "*** Flushing data"
		$PYTHON_EXE $DJANGO_ADMIN flush --settings=$DJANGO_SETTINGS_MODULE > $FIXTURE
	fi

	if [ "$1" = "dump" ]; then
		
		echo "*** Dumping data to $FIXTURE"
		$PYTHON_EXE $DJANGO_ADMIN dumpdata auth.User server --settings=$DJANGO_SETTINGS_MODULE > $FIXTURE
	fi

	if [ "$1" = "migrate" ]; then
		echo "*** Migrating the data into the main database"
		$PYTHON_EXE home/import/migrate.py home/import/datadump
	fi

shift
done
