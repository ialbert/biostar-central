# Setting up the infrastrucuture

Setup is automated via ansible.

Install ansible into the current python environment:

    pip install ansible

Check that ansible works (use your domain instead of test.biostars.org):

    ansible test.biostars.org -u www -m ping all

Alternatively use a host file that contains target server names:

    ansible -i hosts/test.biostars.org -u www -m ping all

## Server configuration

Ensure that you have a generated a public key on your current system. This will be copied over
to the server for public key authentication.

    shh-keygen

Set up the server infrastructure (apt-get, default users, directories)

Manually

    ansible-playbook -i hosts/test.biostars.org server-config.yml

Makefile
     
    make config HOST=hosts/test.biostars.org  

The playbook above will bootstrap a Ubuntu based linux server, a user named `www` and
`nginx` and `postgresql`. You may log into the `www` server via public key authentication.


## Software installation

To install the server and the python dependencies run:

Manually

    ansible-playbook -i hosts/test.biostars.org server-config.yml

Makefile
     
    make install HOST=hosts/test.biostars.org  
    
    
The playbook above will clone the repository into the directory.

    /export/www/biostar-central/
    
At the end of the installation, the playbook will copy the configuration files from

    /export/www/biostar-engine/biostar-central/conf/site/

to

    /export/www/biostar-engine/biostar-central/conf/run/

Link server files manually with:

    ln -sf /export/www/biostar-central/conf/run/site_nginx.conf /etc/nginx/sites-enabled/
    ln -sf /export/www/biostar-central/conf/run/site_supervisor.conf /etc/supervisor/conf.d/

You may now edit and customize the settings file located in `biostar-engine/conf/run/`


Activate and install the conda dependencies:

    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda activate engine
    conda install --file conf/conda_requirements.txt

## Software deployment

To deploy the latest version and restart the servers:

Manually

    ansible-playbook -i hosts/test.biostars.org server-deploy.yml --ask-become-pass

Makefile
     
    make deploy HOST=hosts/test.biostars.org  
    
    
## Migrating from Biostar 1.0

To migrate from an older version of biostar to a server deploying Biostar 2.0

A copy of the older database existing on the local machine can be passed to the playbook. 
This gets copied to the remote server if it does not exist already. 

Manually 
    
    # File path to database we want to migrate
    export LOCAL=/full/path/file.sql.gz
    
    # Migration playbook
	ansible-playbook -i hosts/test.biostars.org server-migrate.yml --ask-become-pass --extra-vars "local_old_db=${LOCAL}"

Makefile 

    make transfer HOST=hosts/test.biostars.org  LOCAL_OLD_DB=/full/path/file.sql.gz
    
