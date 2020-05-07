# Setup Documentation for Bioconductor

Here are steps to running and deploying the forum from scratch.

1. Create a virtual environment and clone most recent version of the forum.


2. Install dependencies. 


3. Run migrations and tests. 


4. Start a local server. 


5. Deploy changes to remote server


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

Activate the `engine` virtual environment.

    $ conda activate engine

Enter the `biostar-central` directory to install dependencies and requirements into the virtual enviorment.

Execute the following to install python requirements: 

    pip install -r conf/requirements.txt      # Install python requirements.
    
    
Add the following conda channels:

    conda config --add channels r
    conda config --add channels conda-forge
    conda config --add channels bioconda

Execute the following to install all anaconda requirements:
    
    conda install --file conf/conda-packages.txt  # Install conda requirements.
    
 After dependencies have been installed, a migration needs to be made to create the database collect static files.
 
 
## 3. Run migrations and tests. 
 
Activate the `engine` virtual environment.

    conda activate engine
    
Migrate the forum app by executing the command:

    python manage.py migrate --settings themes.bioconductor.settings

Collect static files for the forum app by executing the command:

    python manage.py collectstatic --noinput -v 0 --settings themes.bioconductor.settings

There is a `Makefile` command that migrates and collects static files in one shot. 

    make bioconductor init  # Migrate and collect static files. 

A database has now been created and the static files can be found in `biostar-central/export/static/`

To ensure installation and migration was successful, run a test by executing the command: 

    make bioconductor test  # Run tests. 



## Loading demo data and start a local server 

Activate the `engine` virtual environment.

    $ conda activate engine
    
To load sample data into the forum, use the command `make bioconductor startup`:

    make bioconductor startup  # Load sample data
    

Enter the command `make bioconductor serve` to start a local server.

    make bioconductor serve    # Start local server

The site is now available at http://127.0.0.1:8000/. 


You can also load data and start a local server with one command `make bioconductor demo`:

    make bioconductor demo     # Load data and start local server

    
## Default user for local testing site

All users listed in the `ADMINS` attribute of the Django `settings.py` 
module will gain administritave privileges when the site is initialized. The
`DEFAULT_ADMIN_PASSWORD` attribute will be set as the default admin password. 
By default the value for both is:

    admin@localhost

Use this username and password combination to log into the site as an administrator. Change the `DEFAULT_ADMIN_PASSWORD` for public facing installations.


When the site initializes the admin username and password are using the ``ADMINS`` and the ``ADMIN_PASSWORD`` settings in ``biostar/acccounts/settings.py``.

By default both the admin login name and the default admin password are set to

    admin@localhost
   
The Django admin can be found at http://127.0.0.1:8000/accounts/admin/.

## Customize Settings

DO NOT add your custom settings into the public codebase!

The proper practice is to create a separate, independent settings file, then, within that file import **all** default settings. Finally override the fields that you wish to customize in your settings file. For example
create the `my_settings.py` then add into it:

    # Import all default settings.
    from themes.bioconductor.settings import *
    
    # Now override the settings you wish to customize.
    ADMIN_PASSWORD = "foopass"

Apply this settings file with

    python manage.py runserver --settings my_settings.py

Consult the [Django documentation][django] for details.

## Deploying local changes to remote server

Activate the `engine` virtual environment.

    $ conda activate engine

Ensure all of the changes have been committed to the github repo

Enter the following command to deploy local changes to remote server:

    
    make bioconductor deploy  REPO= < github repo url to deploy >       # Deploy to server from a github repo.


### Changing the top banner.

To locally test and edit the top banner, 
you can edit the `themes/bioconductor/templates/banners/top-banners.html`

On the remote site, this top banner is found in `biostar-centra/templates/banners/top-banner.html`

To see the changes in the top-banner.html , restart the server:

    sudo supervisorctl restart forum