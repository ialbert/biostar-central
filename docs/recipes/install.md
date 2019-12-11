
# Getting Started

# Install

The sourcecode can be obtained via::

    git clone https://github.com/ialbert/biostar-central.git

Here are steps to running and deploying the recipes from scratch.

1. Create a virtual environment and clone most recent version.


2. Install dependencies. 


3. Run migrations and tests. 


4. Start a local server. 



## 1. Create environment

Create a virtual environment by first downloading miniconda at https://docs.conda.io/en/latest/miniconda.html. 

After downloading the installation file, run the command ( replace installation_file.sh with your installation file ) : 

    $ bash installation_file.sh      

Once miniconda has been installed, create a virtual enviroment called `engine`.

    $ conda create -n engine python=3.7
    
Start the virtual enviorment by entering the command:

    $ conda activate engine
    
Clone or pull the most recent version of biostars by executing:

      git clone https://github.com/ialbert/biostar-central.git  # Clone a new branch
 
      
## 2. Install dependencies

Activate the `engine` virtual enviorment.

    $ conda activate engine

Enter the `biostar-central` directory to install dependencies and requirements into the virtual enviorment.

Execute the following to install python requirements: 

    pip install -r conf/pip_requirements.txt      # Install python requirements.
    
    
Add the following conda channels:

    conda config --add channels r
    conda config --add channels conda-forge
    conda config --add channels bioconda

Execute the following to install all anaconda requirements:
    
    conda install --file conf/conda_requirements.txt  # Install conda requirements.
    
 After dependencies have been installed, a migration needs to be made to create the database collect static files.
 
 
## 3. Initialize recipes
 
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
    
    make load_recipes
      
    
## 4. Start server 

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

[django]: https://www.djangoproject.com/

# Deploying site

The software follows the recommended practices for developing and deploying [Django web applications][django] .

The [Django documentation][django] contains a wealth of information on the alternative ways to deploy the site on different infrastructure.

Within this setup we recommend the [uwsgi][uwsgi] based deployment.

# Social authentication
