# Setting up the infrastrucuture

The commands below assume an Ubuntu Linux 20.04 LTS distribution. We assume the remote hostname 
and timezones have already been set.

    hostnamectl set-hostname www.foo.com
    timedatectl set-timezone US/Eastern

Setup is automated via [ansible][ansible]. Install ansible into the current python environment:

    pip install ansible

[ansible]: https://www.ansible.com/

## Hosts

The `hosts.ini` file lists the groups of servers that can be targeted. 

Create a copy of this file and add your own hostnames to it. 

The installation commands below will target subsets in hosts file.

## Server Setup

The `server-setup.yml` playbook is designed for a Ubuntu 20.04 LTS based linux server. It will install `nginx`, `postgresql` and other packages, and create the user  `www` that will own the application server install.

Run the playbook with ansible:

    ansible-playbook -i hosts.ini -l test server-setup.yml

The same can be done using `make`:
     
    make setup TARGET=test 

You may need to (manually) restart the server to load some of the updated packages:
    
    reboot now
    
## Software Installation

The ansible playbooks below will perform the following actions:

1. download and install conda, 
2. create a conda enviroment called `engine` prepared to run the biostars software
3. clone the application server and create copies for each configuration file.
4. create copies for the migration and backup scripts.

Run the playbook with ansible:

    ansible-playbook -i hosts.ini -l test server-install.yml

The same can be done using `make`:
     
    make install TARGET=test   
    
The playbook above will clone the repository into the directory.

    /export/www/biostar-central/
    
At the end of the installation, the playbook will copy the configuration files from

    /export/www/biostar-engine/biostar-central/conf/site/

to 

    /export/www/biostar-engine/biostar-central/conf/run/

You will need to link the configuration files from command line:

    ln -sf /export/www/biostar-central/conf/run/site_nginx.conf /etc/nginx/sites-enabled/
    ln -sf /export/www/biostar-central/conf/run/site_supervisor.conf /etc/supervisor/conf.d/

You may now edit and customize the settings file located in `biostar-engine/conf/run/`

By default only the immediate software requirements are installed. 

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
    
