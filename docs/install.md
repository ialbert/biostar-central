## Documentation for Bioinformatics Recipes

### What is the default admin login?

When the site initializes the admin username and password are using the `ADMINS` and the `ADMIN_PASSWORD` settings in `biostar/acccounts/settings.py`.
 
 By default both the admin login name and the default admin password are set to

    admin@localhost

**Note**: These settings must be changed on a publicly accessible site!

### How to access the Django Admin interface?

* http://127.0.0.1:8000/accounts/admin/

### How to customize the settings?

DO NOT add your custom settings into the public codebase!

The proper practice is to create a separate, independent settings file, then, within that file import **all** default settings. Finally override the fields that you wish to customize in your settings file. For example
create the `my_settings.py` then add into it:

    # Import all default settings.
    from biostar.recipes.settings import *

    # Now override the settings you wish to customize.
    ADMIN_PASSWORD = "foopass"

Apply this settings file with

    python manage.py runserver --settings my_settings.py