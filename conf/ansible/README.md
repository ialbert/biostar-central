# System setup

Setup is automated via [ansible][ansible]. Install ansible into the current python environment:
    
    # Using pip
    pip install ansible
    
    # Using apt-get
    sudo apt install ansible 
    
The commands below assume the following:

1. An Ubuntu Linux 20.04 LTS distribution. 
2. You are able to log into the `root` user of the host (public-key authentication preferred)
2. The remote hostname and timezones have already been set.

        hostnamectl set-hostname www.foo.com
        timedatectl set-timezone US/Eastern

[ansible]: https://www.ansible.com/

## Hosts

The `hosts.ini` file lists the groups of servers that can be targeted. For example:

    [test]
    test.biostars.org
    
    [recipes]
    www.bioinformatics.recipes

The installation commands below will target subsets in hosts file.

## Quick start

Run:

    make setup install deploy TARGET=test 

Once completed a default site will be installed and deployed via Nginx and Postgresql.

Read on for details on what takes place during each step.

## Server setup

The `server-setup.yml` playbook is designed for a Ubuntu 20.04 LTS based linux server. It will install `nginx`, `postgresql` and other packages, creates the user  `www` that will own the application server install.

    make setup TARGET=test 
    
or you can run the playbook with ansible:

    ansible-playbook -i hosts.ini -l test server-setup.yml

You may need to manually restart the server to apply some of the updates:
    
    reboot now
    
## Software install

The ansible playbooks below will perform the following actions as user `www`:

1. download and install `conda`, 
1. create a conda enviroment called `engine` that can run the biostars software
1. clone the source code for the application server 
1. create local copies for the configurations (`conf/run` folder)
1. create data migration and backup scripts 

To perform the install run:
     
    make install TARGET=test   

Or you may run the ansible playbook directly:

    ansible-playbook -i hosts.ini -l test server-install.yml
    
The playbook above will clone the repository into the directory.

    /export/www/biostar-central/
    
Edit the settings file located in `conf/run/site_settings.py` and change:

* `SECRET_KEY`
* `SITE_DOMAIN`

Add more settings as needed then restart the servers:

    make restart TARGET=test
    
## Database settings

By default the postgresql database will be accessible only from the local 
machine with `www` user having database creation roles.

## Sending mail

The default site set up with ansible will send email via SMTP. To set up your system run:

    aptitude install mailutils

## Software deployment

To deploy the latest version and restart the servers:
 
    make deploy TARGET=test  
    
    make deploy TARGET=demo REPO=https://github.com/Bioconductor/SupportUpgrade/tree/master/conf USER=ubuntu
        
or via the playbook:

    ansible-playbook -i hosts/test.biostars.org server-deploy.yml --ask-become-pass

## Restart remote 

    make restart TARGET=test
   

## Deploying from different branches

To deploy from different branches, log into remote server:

    git checkout < branch >
    

Then locally run `make deploy TARGET=<target>` with the target set to the remote server.

This will apply all changes found in a branch.
 
## Migrating from Biostar 1.0 (TODO)

To migrate from an older version of biostar to a server deploying Biostar 2.0

A copy of the older database existing on the local machine can be passed to the playbook. 
This gets copied to the remote server if it does not exist already. 

Manually 
    
    # File path to database we want to migrate
    export LOCAL=/full/path/file.sql.gz
    
    # Migration playbook
	ansible-playbook -i test.biostars.org server-migrate.yml --extra-vars "local_old=${LOCAL}"

Makefile 

    make transfer TARGET=test LOCAL_OLD_DB=/full/path/file.sql.gz
 
Deploy and transfer using different user,host,repo, and data dump.

    make setup install deploy transfer TARGET=demo USER=ubuntu REPO=https://github.com/Natay/biostar-central.git LOCAL_OLD_DB=~/tmp/data-dump.gz