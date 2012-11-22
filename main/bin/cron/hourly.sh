#!/bin/bash
echo "*** hourly.sh: `date`"
source /home/biostar/www/production/live/biostar.env

set -ue

$PYTHON_EXE -m main.server.search

