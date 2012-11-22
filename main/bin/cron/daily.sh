#!/bin/bash
echo "*** daily: `date`"
source /home/biostar/www/production/live/biostar.env

set -ue

$PYTHON_EXE -m main.bin.sitemap
$PYTHON_EXE -m main.bin.planet --download
$PYTHON_EXE -m main.bin.patch --reduce_notes -n 6
