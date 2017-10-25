
# Make sure Main results page doesn't fail out.
touch {{settings.index}}

cd ../../../..

python manage.py data --copy --fname1 {{copy_from.path}} --fname2 {{copy_to.path}}