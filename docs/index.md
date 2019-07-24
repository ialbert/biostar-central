# Biostar Central Documentation

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

* [recipe-concepts.md](recipe-concepts.md)
* [recipe-commands.md](recipe-commands.md)
* [recipe-deploy.md](recipe-deploy.md)
* [recipe-api.md](recipe-api.md)

## How to set up the `forum` app

(TODO)

## Django configuration

The software follows the recommended practices for developing and deploying [Django web applications][django] .

The [Django documentation][django] contains a wealth of information on the alternative ways to deploy the site on different infrastructure.

## Infrastructure configuration

To run bioinformatics oriented software via the recipes additional configuration of your environment may be necessary. For example we use:

    conda config --add channels r
    conda config --add channels conda-forge
    conda config --add channels bioconda

    # Install the conda requirements.
    conda install --file conf/conda_requirements.txt
