# Setup is automated via ansible

Install ansible into the current python environment:

    pip install ansible

Check that ansible works (use your domain instead of test.biostars.org) :

    ansible test.biostars.org -u www -m ping all

or using a host file that contains target server names:

    ansible -i hosts/test.biostars.org -u www -m ping all

## Set up the infrastructure

Ensure that you have a generated a public key on your current system. This will be copied over
to the server for public key authentication.

    shh-keygen

Set up the server infrastructure (apt-get, default users, directories)

    ansible-playbook -i hosts/test.biostars.org server-config.yml

Install the software:

    ansible-playbook -i hosts/test.biostars.org software-install.yml

By default conda is not added to the `PATH`, is that behavior is desired one needs to run:

    ~/miniconda3/bin/conda init

Link the nginx and supervisor configurations:

    ln -sf /export/sites/biostar-engine/conf/site/site_nginx.conf /etc/nginx/sites-enabled/
    ln -sf /export/sites/biostar-engine/conf/site/site_supervisor.conf /etc/supervisor/conf.d/ 
        
Initialize the certificates

    sudo certbot --nginx

You should test your configuration at:

* https://www.ssllabs.com/ssltest/analyze.html?d=bioinformatics.recipes
* https://www.ssllabs.com/ssltest/analyze.html?d=data.bioinformatics.recipes
* https://www.ssllabs.com/ssltest/analyze.html?d=www.bioinformatics.recipes

## Deployment

The following will pull the new content and restart the servers:

    ansible-playbook -i hosts server_deploy.yml --ask-become-pass


To restart servers alone:

    ansible-playbook -i hosts server_deploy.yml --ask-become-pass --extra-vars "restart=True"


To install dependencies:

    ansible-playbook -i hosts server_deploy.yml --ask-become-pass --extra-vars "install=True restart=True"

To reset the site:

    ansible-playbook -i hosts server_deploy.yml --ask-become-pass --extra-vars "reset=True"

