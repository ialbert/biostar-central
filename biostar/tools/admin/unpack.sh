
touch {{settings.index}}
cd ../../../..
FNAME={{unpack.path}}

python manage.py data --unpack --fname1=${FNAME}

