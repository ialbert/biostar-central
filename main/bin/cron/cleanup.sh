#!/bin/bash
echo "*** cleanup.sh: `date`"
source /home/biostar/www/production/live/biostar.env

set -ue

$PYTHON_EXE $DJANGO_ADMIN cleanup
#$PYTHON_EXE $DJANGO_ADMIN openid_cleanup

