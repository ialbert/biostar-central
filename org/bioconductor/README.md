# Bioconductor Support Site

How to customize Biostar

## Custom settings

All custom settings should import from the main settings:

    from biostar.settings.base import *
    
Within this new file users need to only override only the settings that are to be changed.

## Customizing templates

Templates that are to be modified need to be copied to a new directory. Example from:

    biostar/server/templates 
    
directory to:

    org/bioconductor/templates

so that it retains the same relative path.  If the template was located
in a subdirectory of the `templates` directory this same
subdirectory would need to be created in the new location as well.

## Template paths

To activate the Bioconductor templates users need to ensure that

    org/bioconductor/templates

directory is listed ahead of the other directories in `TEMPLATE_DIRS` variable.
For example in the custom settings module:
    
    TEMPLATE_DIRS = [THEME_PATH] + TEMPLATE_DIRS

## Enabling the new settings

Run biostar in an environment where the new file is listed as a python module name 
(not file name) or the `DJANGO_SETTINGS_MODULE` environment variable

    # The file is called foo.py
    export DJANGO_SETTINGS_MODULE=foo

The settings file must be located in a directory that Python can import. This can 
be achieved by adding the directory that stores the new settings file to the Python path:
    
    # The full path is /home/biostar/bar/foo.py
    export PYTHONPATH=/home/biostar/bar:$PYTHONPATH
    export DJANGO_SETTINGS_MODULE=foo

## Environment variables

Biostar may be customized via a number of environment variables that otherwise 
have default values. Export an environment variable with the given name to
override a setting. Default values:

    DJANGO_SETTINGS_MODULE = biostar.settings.base
    BIOSTAR_ADMIN_NAME = Biostar Admin
    BIOSTAR_ADMIN_EMAIL = admin@lvh.me
    SECRET_KEY = admin@lvh.me
    EMAIL_HOST = localhost
    EMAIL_PORT = 25
    EMAIL_HOST_USER = postmaster
    EMAIL_HOST_PASSWORD = password

These settings are for convenience only. The custom settings file that a user
creates may choose to directly specify any of these settings.