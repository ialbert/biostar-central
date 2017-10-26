
echo Redirecting
ID={{target.value}}

(cd ../../../.. && python manage.py admintasks --unpack --target=${ID})

