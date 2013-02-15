#!/bin/sh

if [ -z "$BIOSTAR_HOME" ]; then
    echo "(!) environment variables not set"
    echo "Try: source conf/default.env"
    exit 1
fi

GUNICORN=/export/2.7/bin/gunicorn
ROOT=/Users/ialbert/app/biostar-central
PID=/Users/ialbert/app/biostar-central/apache/logs/gunicorn.pid
CONF=$ROOT/apache/gunicorn.conf.py

APP=conf.test_deploy_wsgi:application

if [ -f $PID ]; then rm $PID; fi

cd $ROOT
echo "*** starting gunicorn with config CONF=$CONF, PID=$PID, APP=$APP"
exec $GUNICORN -c $CONF --pid=$PID $APP
