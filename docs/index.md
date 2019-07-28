# Biostar Central Documentation

## How to use the `recipes`?

* [Recipe concepts](recipe-concepts.md)
* [recipe-commands.md](recipe-commands.md)
* [recipe-deploy.md](recipe-deploy.md)
* [recipe-api.md](recipe-api.md)

## How to customize the `forum`?

(TODO)

## What is the default admin password?

The admin username and password are set via the `ADMINS` and the `ADMIN_PASSWORD` settings in `biostar/accounts/settings.py`. By default both the admin login name and the default admin password are set to

    admin@localhost

**Note**: These settings must be changed on a publicly facing site!

## How to customize the settings?

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

## How to configure Django?

The software follows the recommended practices for developing and deploying [Django web applications][django] .

The [Django documentation][django] contains a wealth of information on the alternative ways to deploy the site on different infrastructure.

## How to configure the computer?

To run bioinformatics oriented software via the recipes additional configuration of your environment may be necessary. For example we use:

    conda config --add channels r
    conda config --add channels conda-forge
    conda config --add channels bioconda

    # Install the conda requirements.
    conda install --file conf/conda_requirements.txt
