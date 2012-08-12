#!/bin/bash

if [ -z "$BIOSTAR_HOME" ]; then
    echo "(!) environment variables not set"
    echo "Try: source conf/default.env"
    exit 1
fi

set -ue

# verbosity level for commands 0=minimal, 2=maximal
VERBOSITY=1

# the migration path and limit
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
        echo "*** PG_USER=$PG_USER"

        echo "--- migration"
        echo "*** MIGRATE_PATH=$MIGRATE_PATH"
        echo "*** MIGRATE_LIMIT=$MIGRATE_LIMIT"
        echo "*** JSON_FIXTURE=$JSON_FIXTURE"
        
        echo "--- environment "
        echo "*** DJANGO_ADMIN=$DJANGO_ADMIN"
        echo "*** PYTHON_EXE=$PYTHON_EXE"
        echo "*** BIOSTAR_HOME=$BIOSTAR_HOME"
        echo "*** BIOSTAR_HOSTNAME=$BIOSTAR_HOSTNAME"
        echo "*** VERBOSITY=$VERBOSITY"    
    fi
    
    if [ "$1" = "delete" ]; then
        # deletes the sqlite database
        if [ -n "$SQLITE_DBNAME" ]; then
            echo "*** deleting SQLITE_DBNAME=$SQLITE_DBNAME"
            rm -f $SQLITE_DBNAME
        else
            echo "(!) SQLITE_DBNAME variable is empty"
        fi
    fi
 
    if [ "$1" = "planet" ]; then
        # initializes the planet
        echo '*** initializes the planet'
        $PYTHON_EXE -m main.scripts.planet --init 30 --download --update 1
    fi

    if [ "$1" = "flush" ]; then
        echo "*** flushing the database"
        $PYTHON_EXE $DJANGO_ADMIN flush --noinput --settings=$DJANGO_SETTINGS_MODULE
    fi

    if [ "$1" = "init" ]; then
        echo "*** initializing server on $BIOSTAR_HOSTNAME"
        $PYTHON_EXE $DJANGO_ADMIN syncdb -v $VERBOSITY --noinput --settings=$DJANGO_SETTINGS_MODULE
        $PYTHON_EXE $DJANGO_ADMIN migrate main.server --settings=$DJANGO_SETTINGS_MODULE
         
        echo "*** collecting static files"
        $PYTHON_EXE $DJANGO_ADMIN collectstatic -v $VERBOSITY --noinput --settings=$DJANGO_SETTINGS_MODULE
    fi

    if [ "$1" = "import" ]; then
        echo "*** importing data from $JSON_FIXTURE"
        $PYTHON_EXE $DJANGO_ADMIN loaddata $JSON_FIXTURE --settings=$DJANGO_SETTINGS_MODULE
    fi
    
    if [ "$1" = "dump" ]; then        
        echo "*** dumping data to $JSON_FIXTURE"
        $PYTHON_EXE $DJANGO_ADMIN dumpdata auth.User server --settings=$DJANGO_SETTINGS_MODULE | gzip > $JSON_FIXTURE
    fi

    if [ "$1" = "pgcreate" ]; then
        # creates the PG database
        echo "*** creating postgresql database PG_DBNAME=$PG_DBNAME"
        createdb $PG_DBNAME
    fi
    
    if [ "$1" = "pgdrop" ]; then
        # drops the PG database
        echo "*** dropping postgresql database PG_DBNAME=$PG_DBNAME"
        dropdb $PG_DBNAME
    fi
    
    if [ "$1" = "pgreset" ]; then
        # resets the postgresql database, removes all tables
        echo '*** create drop table commands'
        $PYTHON_EXE $DJANGO_ADMIN sqlclear server django_openid_auth sites sessions admin auth contenttypes > import/sqlclear.sql --settings=$DJANGO_SETTINGS_MODULE
        echo '*** postgresql reset'
        psql -U $PG_USER $PG_DBNAME < import/sqlclear.sql
    fi
    
    if [ "$1" = "pgdump" ]; then
        # dumps a postgres database to a file
        echo "*** dumping database $PG_DBNAME"
        pg_dump -O -x $PG_DBNAME 
    fi
    
     if [ "$1" = "pgimport" ]; then
        # restores a postgresl database from a file
        echo "*** restoring database $PG_DBNAME"
        psql -U $PG_USER $PG_DBNAME < $2
    fi
    
    if [ "$1" = "test" ]; then
        echo "*** running the tests"
        #$PYTHON_EXE $DJANGO_ADMIN test server --settings=$DJANGO_SETTINGS_MODULE --failfast
        $PYTHON_EXE $DJANGO_ADMIN test server --settings=$DJANGO_SETTINGS_MODULE --failfast
    fi

    if [ "$1" = "migrate" ]; then
        echo "*** MIGRATE_PATH=$MIGRATE_PATH"
        echo "*** MIGRATE_LIMIT=$MIGRATE_LIMIT"
        $PYTHON_EXE -m main.scripts.migrate --path $MIGRATE_PATH --limit $MIGRATE_LIMIT
    fi

    if [ "$1" = "index" ]; then        
        echo "*** indexing all post content"
        $PYTHON_EXE -m main.server.search
    fi

    if [ "$1" = "run" ]; then
        echo "*** running the webserver on $BIOSTAR_HOSTNAME"
        $PYTHON_EXE $DJANGO_ADMIN runserver $BIOSTAR_HOSTNAME --settings=$DJANGO_SETTINGS_MODULE
    fi
    
    if [ "$1" = "selenium" ]; then
        # needs to reindex to be most up to date
        $PYTHON_EXE -m main.server.search
        
        echo "*** running selenium on $SELENIUM_TEST_URL"
        $PYTHON_EXE main/server/tests/selenium_tests.py $SELENIUM_TEST_URL
    fi
    
    if [ "$1" = "deploy" ]; then
        echo "*** deploying biostar the the remote server"
        
        # read off the deployment variables
        source conf/deploy.env
        
        # remove the index
        rm -rf main/db/index/*

        # migrate the entire datadump and initialize the planet
        time ./biostar.sh pgdrop init migrate planet index
        
        # dump the to an SQL fixture
        echo "--- generating $SQL_FIXTURE"
        time ./biostar.sh pgdump
        
        # upload the fixture to the remote
        echo "--- uploading $SQL_FIXTURE to $REMOTE_IMPORT"
        rsync -azv $SQL_FIXTURE $REMOTE_IMPORT
        
        # synchronize the remote index
        echo "--- uploading the indices"
        rsync -azv main/db/index/* $REMOTE_INDEX

        echo ""
        echo "--- uploaded $SQL_FIXTURE to the remote location $REMOTE_IMPORT"
        echo ""
    fi
       
shift
done

