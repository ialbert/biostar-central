#!/bin/bash
echo "*** sixhours.sh: `date`"
source /home/biostar/www/production/live/biostar.env

set -ue

$PYTHON_EXE -m main.bin.planet --update 1
$PYTHON_EXE -m main.bin.patch --blog_cleanup
