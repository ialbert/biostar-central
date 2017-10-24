
touch {{settings.index}}
cd ../../../..
python manage.py test_email --to {{email.value}}