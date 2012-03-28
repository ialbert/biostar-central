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

export SQLITE_DBNAME=${SQLITE_DBNAME:-"$BIOSTAR_HOME/db/biostar.db"}

# the migration path and limit
export MIGRATE_PATH=${MIGRATE_PATH:-"import/se0"}
export MIGRATE_LIMIT=${MIGRATE_LIMIT:-"100"}

# the JSON fixture to dump/load data from
export JSON_FIXTURE=${JSON_FIXTURE:-"import/biostar-dump.json.gz"}

# the SQL fixture to dump/load postgresql SQL 
export SQL_FIXTURE=${SQL_FIXTURE:-"import/biostar-dump.sql"}

# the postgresql database name
export PG_DBNAME=${PG_DBNAME:-"biostar-test-database"}

# the DJANGO_SETTINGS_MODULE needs to be in the python import path
export PYTHONPATH=$PYTHONPATH:$BIOSTAR_HOME   

# add the library files to the pythonpath
export PYTHONPATH=$PYTHONPATH:libs/:libs/libraries.zip

# setting up the python
export PYTHON_EXE=${PYTHON_EXE:-"python"}
export DJANGO_ADMIN=manage.py

echo "--- main settings"
echo "*** BIOSTAR_HOME=$BIOSTAR_HOME"
echo "*** BIOSTAR_HOSTNAME=$BIOSTAR_HOSTNAME"
echo "*** DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE"
#echo "*** PYTHONPATH=$PYTHONPATH"

if [ $# == 0 ]; then
	echo ''
	echo 'Usage:'
	echo '  $ run.sh <command>'
	echo ''
	echo 'Multiple commands may be used in the same line:'
	echo '  $ run.sh init import run'
	echo ''
	echo 'Commands:'
	echo '  init     - initializes the database'
	echo '  import   - imports a data fixture'
    echo '  dump     - dumps the current database as a data fixture'
	echo '  delete   - removes the sqlite database (sqlite specific)'
	echo '  run      - runs server'
	echo '  migrate  - parses a StackExchange XML dump to a data fixture'
	echo '  test     - runs all tests'
    echo '  env      - shows all customizable environment variables'
    echo ''
    echo 'Use environment variables to customize the behavior. See docs.'
fi

while (( "$#" )); do

	if [ "$1" = "env" ]; then
        echo "--- databases"
        echo "*** SQLITE_DBNAME=$SQLITE_DBNAME"
        echo "*** PG_DBNAME=$PG_DBNAME"
        
        echo "--- migration"
        echo "*** MIGRATE_PATH=$MIGRATE_PATH"
        echo "*** MIGRATE_LIMIT=$MIGRATE_LIMIT"
        echo "*** JSON_FIXTURE=$JSON_FIXTURE"
        echo "*** SQL_FIXTURE=$SQL_FIXTURE"
        
        echo "--- environment "
        echo "*** DJANGO_ADMIN=$DJANGO_ADMIN"
        echo "*** PYTHON_EXE=$PYTHON_EXE"
        echo "*** BIOSTAR_HOME=$BIOSTAR_HOME"
        echo "*** BIOSTAR_HOSTNAME=$BIOSTAR_HOSTNAME"
        echo "*** VERBOSITY=$VERBOSITY"    
	fi
    
	if [ "$1" = "delete" ]; then
        # deletes the sqlite database
		echo "*** deleting sqlite"
		rm -f $SQLITE_DBNAME
	fi

    if [ "$1" = "pgdrop" ]; then
        # drops the PG datanase
        echo '*** dropping postgresql'
        dropdb $PG_DBNAME
        echo '*** creating postgresql'
        createdb $PG_DBNAME
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
		echo "*** importing data from $JSON_FIXTURE"
		$PYTHON_EXE $DJANGO_ADMIN loaddata $JSON_FIXTURE --settings=$DJANGO_SETTINGS_MODULE
	fi

    if [ "$1" = "pgimport" ]; then
        # restores a postgresl database from a file
        echo "*** restoring database $PG_DBNAME from $SQL_FIXTURE"
        psql $PG_DBNAME < $SQL_FIXTURE
	fi
    
	if [ "$1" = "test" ]; then
		echo "*** running the tests"
		#$PYTHON_EXE $DJANGO_ADMIN test server --settings=$DJANGO_SETTINGS_MODULE --failfast
		$PYTHON_EXE $DJANGO_ADMIN test server --settings=$DJANGO_SETTINGS_MODULE --failfast
	fi

    if [ "$1" = "pgdump" ]; then
        # dumps a postgres database to a file
        echo "*** dumping database $PG_DBNAME to $SQL_FIXTURE"
        pg_dump $PG_DBNAME > $SQL_FIXTURE
        wc -l $SQL_FIXTURE
	fi
    
	if [ "$1" = "dump" ]; then		
		echo "*** dumping data to $JSON_FIXTURE"
		$PYTHON_EXE $DJANGO_ADMIN dumpdata auth.User server --settings=$DJANGO_SETTINGS_MODULE | gzip > $JSON_FIXTURE
	fi

	if [ "$1" = "migrate" ]; then
		echo "*** migrating data to a $JSON_FIXTURE"
		echo "*** MIGRATE_PATH=$MIGRATE_PATH"
		echo "*** MIGRATE_LIMIT=$MIGRATE_LIMIT"
		echo "*** JSON_FIXTURE=$JSON_FIXTURE"
		$PYTHON_EXE -m main.scripts.migrate --path $MIGRATE_PATH --limit $MIGRATE_LIMIT -o $JSON_FIXTURE.temp
        cat $JSON_FIXTURE.temp | gzip > $JSON_FIXTURE
        rm $JSON_FIXTURE.temp
        echo "*** migrated data to a $JSON_FIXTURE"
	fi

	if [ "$1" = "index" ]; then		
		echo "*** indexing all post content"
		$PYTHON_EXE -m main.server.search
	fi

    if [ "$1" = "run" ]; then
		echo "*** running the webserver on $BIOSTAR_HOSTNAME"
		$PYTHON_EXE $DJANGO_ADMIN runserver $BIOSTAR_HOSTNAME --settings=$DJANGO_SETTINGS_MODULE
	fi
    
shift
done

