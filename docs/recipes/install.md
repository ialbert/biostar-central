# Installation

The code in **Biostar Recipes**  requires [Python 3.6][python] or above.

Our installation instructions rely on [conda][conda] though other alternatives for managing python environments are equally viable.

    # Create a virtual environment.
    conda create -y --name engine python=3.6

    # Activate the python environment.
    conda activate engine

    # Clone the source server code and the recipe code.
    git clone https://github.com/ialbert/biostar-central.git

    # Switch to the biostar-central directory.
    cd biostar-central

    # Install server dependencies.
    pip install -r conf/requirements.txt

The installation is now complete. All server management commands are run through `make` by running one or more `make` tasks.
For example to test the `recipes` app run:

    make recipes test

## Running a Demo

To run the demonstration version of the `recipes` app execute:

    make recipes demo

Visit <http://127.0.0.1:8000/> to view the site.

##  Initialize Recipes
 
Activate the `engine` virtual enviorment.

    conda activate engine
    
Migrate the recipes app by executing the command:

    python manage.py migrate --settings biostar.recipes.settings

Collect static files for the recipes app by executing the command:

    python manage.py collectstatic --noinput -v 0 --settings biostar.recipes.settings

There is a `Makefile` command that migrates and collects static files in one shot. 

    make recipes init  # Migrate and collect static files. 

A database has now been created and the static files can be found in `biostar-central/export/static/`

To ensure installation and migration was successful, run a test by executing the command: 

    make recipes test  # Run tests. 
    
  
To populate the database with random data run:
    
    make recipes startup
      
    
## Start Server 

Activate the `engine` virtual enviorment:

    $ conda activate engine
    
Start a local server:

    make recipes serve    # Start local server

The site is now available at http://127.0.0.1:8000/. 
 
When the site initializes the admin username and password are using the ``ADMINS`` and the ``ADMIN_PASSWORD`` settings in ``biostar/acccounts/settings.py``.

By default both the admin login name and the default admin password are set to

    admin@localhost
   
The Django admin can be found at http://127.0.0.1:8000/accounts/admin/.

# Customize Settings

DO NOT add your custom settings into the public codebase!

The proper practice is to create a separate, independent settings file, then, within that file import **all** default settings. Finally override the fields that you wish to customize in your settings file. For example
create the `my_settings.py` then add into it:

    # Import all default settings.
    from biostar.recipes.settings import *

    # Now override the settings you wish to customize.
    ADMIN_PASSWORD = "foopass"

Apply this settings file with

    python manage.py runserver --settings my_settings.py

Consult the [Django documentation][django] for details.


## Directory Structure 

Each project has a physical directory associated on the system located on the system. 


1. Projects directory
    - Each project has a directory with the data associated. 
2. Results directory
   - Location where the results of a recipe run are stored.
3. Table of contents directory
    - Contains table of content files for every data.

These directories all found in the media directory found in the `settings.py` under `MEDIA_ROOT`. The general structure is:

    media/
        projects/
           ...
        jobs/
           ...
        tocs/
            ... 

[django]: https://www.djangoproject.com/

# Deploying Site

The software follows the recommended practices for developing and deploying [Django web applications][django] .

The [Django documentation][django] contains a wealth of information on the alternative ways to deploy the site on different infrastructure.

Within this setup we recommend the [uwsgi][uwsgi] based deployment.

[python]: https://www.python.org/
[django]: https://www.djangoproject.com/
[biostars]: https://www.biostars.org
[recipes]: https://www.bioinformatics.recipes
[handbook]: https://www.biostarhandbook.com
[conda]: https://conda.io/docs/
