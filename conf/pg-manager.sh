!/bin/bash
set -ue

# verbosity level for commands 0=minimal, 2=maximal
VERBOSITY=1

# should be used from the main BIOSTAR home directory

SQL_DUMP=${SQL_DUMP:-"import/datadump.sql"}
DB_NAME=${DB_NAME:-"biostar-test-database"}
echo ""
echo "Settings:"
echo "*** SQL_DUMP=$SQL_DUMP"
echo "*** DB_NAME=$DB_NAME"
echo "*** DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE"

if [ $# == 0 ]; then
	echo ''
	echo Usage:
	echo '  $ pg-manager.sh <command>'
	echo ''
	echo 'Multiple commands may be used in the same line:'
	echo '  $ run.sh dump restore'
	echo ''
	echo 'Commands:'
	echo '  dump     - dumps a database into a file'
	echo '  restore  - restores database from a file'
fi

while (( "$#" )); do

    if [ "$1" = "dump" ]; then
        # dumps the database to a file
        echo "*** dumping database $DB_NAME to $SQL_DUMP"
        pg_dump $DB_NAME > $SQL_DUMP
        wc -l $SQL_DUMP
	fi
    
     if [ "$1" = "restore" ]; then
        # restores data from a file
        echo "*** restoring database $DB_NAME from $SQL_DUMP"
        echo '*** dropping postgresql'
        dropdb $DB_NAME
        echo '*** creating postgresql'
        createdb $DB_NAME
        psql $DB_NAME < $SQL_DUMP
	fi
    
shift
done