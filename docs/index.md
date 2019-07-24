
# Biostar Central documentation

### Important note

DO NOT add your custom settings into the public codebase!

The proper practice is to create a separate, independent settings file, then, within that file import **all** default settings. Finally override the fields that you wish to customize in your settings filr. For example
create the `foo_settings.py` then add into it:

    # Import all default settings.
    from biostar.recipes.settings import *

    # Now override the settings you wish to customize.
    ADMIN_PASSWORD = "foopass"

Apply this settings file with

    python manage.py runserver --settings foo_settings.py

Consult the The [Django documentation][django] for details.

## How to set up the `recipes` app

* [recipe-commands.md](recipe-commands.md)


## How to set up the `forum` app

