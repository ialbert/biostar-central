# Install

The sourcecode can be obtained via::

    git clone https://github.com/ialbert/biostar-central.git

Here are steps to running and deploying the forum from scratch.


1. Create a virtual environment and clone most recent version of the forum.


2. Install dependencies. 


3. Run migrations and tests. 


4. Start a local server. 


## 1. Create a virtual environment and clone the repo.

Create a virtual environment by first downloading miniconda at https://docs.conda.io/en/latest/miniconda.html. 

After downloading the installation file, run the command ( replace installation_file.sh with your installation file ) : 

    $ bash installation_file.sh      

Once miniconda has been installed, create a virtual enviroment called `engine`.

    $ conda create -n engine python=3.7
    
Start the virtual enviorment by entering the command:

    $ conda activate engine
    
Clone or pull the most recent version of the forum by executing:

      git clone https://github.com/ialbert/biostar-central.git  # Clone a new branch
      
      git pull https://github.com/ialbert/biostar-central.git   # Pull into an exisiting 
      
      
## 2. Install dependencies. 

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
 
 
 ## 3. Run migrations and tests. 
 
Activate the `engine` virtual enviorment.

    conda activate engine
    
Migrate the forum app by executing the command:

    python manage.py migrate --settings biostar.forum.settings

Collect static files for the forum app by executing the command:

    python manage.py collectstatic --noinput -v 0 --settings biostar.forum.settings

There is a `Makefile` command that migrates and collects static files in one shot. 

    make forum init  # Migrate and collect static files. 

A database has now been created and the static files can be found in `biostar-central/export/static/`

To ensure installation and migration was successful, run a test by executing the command: 

    make forum test  # Run tests. 
    
    
## 4. Start a local server 

Activate the `engine` virtual enviorment.

    $ conda activate engine
    
Enter the command `make forum serve` to start a local server.

    make forum serve    # Start local server

The site is now available at http://127.0.0.1:8000/. 
 


## Default Login

When the site initializes the admin username and password are using the ``ADMINS`` and the ``ADMIN_PASSWORD`` settings in ``biostar/acccounts/settings.py``.

 By default both the admin login name and the default admin password are set to

.. code-block::

   admin@localhost
   
The sourcecode can be obtained via::

	git clone https://github.com/ialbert/biostar-central.git


### How to access the Django Admin interface?

* http://127.0.0.1:8000/accounts/admin/

## How to customize the settings?

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

## How do I deploy the site?

The software follows the recommended practices for developing and deploying [Django web applications][django] .

The [Django documentation][django] contains a wealth of information on the alternative ways to deploy the site on different infrastructure.

Within this setup we recommend the [uwsgi][uwsgi] based deployment.

## Social authentication


## Migrating from Biostar 1.X
