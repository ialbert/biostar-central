## Customize

**Note:** when customizing any part of Biostar users do not typically need to modify the original
files. Typically one would create a new folders, files or settings and instruct Biostar to load 
these **in addition** to the original file. When looking for any file 
there is a priority of the search order. Any file with the same name that is found before
the original file will override it. Using the site this way will allow changes from the 
main trunk to be more easily incorporated.

### Customizing settings 

[django]: https://www.djangoproject.com/

Biostar is a [Django][django] based application and follows the Django conventions for loading
settings. Upon starting up the value of the `DJANGO_SETTINGS_MODULE` environment variable determines
which settings module will be imported.

To more easily change these settings the default Biostar configuration will also read more
environment variables.






