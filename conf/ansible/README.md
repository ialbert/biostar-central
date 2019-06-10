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

    ansible-playbook -i hosts/test.biostars.org server-config.yml

The playbook above will bootstrap a Ubuntu based linux server, a user named `www` and
`nginx` and `postgresql`. You may log into the `www` server via public key authentication.

## Software installation

To install the server and the python dependencies run:

    ansible-playbook -i hosts/test.biostars.org software-install.yml

The playbook above will clone the repository into the directory.

    /export/www/biostar-engine/

At the end of the installation, the playbook will copy the configuration files from

    /export/www/biostar-engine/biostar-engine/conf/site/

to

    /export/www/biostar-engine/biostar-engine/conf/run/

Link server files manually with:

    ln -sf /export/www/biostar-engine/conf/run/site_nginx.conf /etc/nginx/sites-enabled/
    ln -sf /export/www/biostar-engine/conf/run/site_supervisor.conf /etc/supervisor/conf.d/

You may now edit and customize the settings file located in `biostar-engine/conf/run/`

By default conda will not be added to the `PATH` variable of the `www` user. If that behavior is desired one needs to run:

    ~/miniconda3/bin/conda init

Activate and install the conda dependencies:

    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda activate engine
    conda install --file conf/conda_requirements.txt

## Software deployment

To deploy the latest version and restart the servers:

    ansible-playbook -i hosts/test.biostars.org server-deploy.yml --ask-become-pass

To restart servers alone:

    ansible-playbook -i hosts/test.biostars.org server-deploy.yml --ask-become-pass --extra-vars "restart=True"


To install dependencies:

    ansible-playbook -i hosts/test.biostars.org server-deploy.yml --ask-become-pass --extra-vars "install=True restart=True"

To reset the site:

    ansible-playbook -i hosts/test.biostars.org server-deploy.yml --ask-become-pass --extra-vars "reset=True"

