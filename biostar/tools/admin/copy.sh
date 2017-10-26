
SRC={{source.value}}
TRG={{target.value}}

(cd ../../../.. && python manage.py admintasks --copy --source=${SRC} --target=${TRG})