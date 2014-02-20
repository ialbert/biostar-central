"""
Command run when intializing an Ubuntu based linux distro from scratch
"""

from fabric.api import run, cd
from fabric.context_managers import prefix
from fabric.api import *
from getpass import getpass
from sites import *

def usegalaxy(user="www"):
    "Sets the environment for the biostar galaxy"
    env.hosts = ['biostar.usegalaxy.org']
    env.user = user

def add_ssh_key():
    "Appends the current SSH pub key to the remote authorized keys"
    local("scp ~/.ssh/id_rsa.pub %(user)s@%(host)s:~/" % env)
    run("cat ~/id_rsa.pub >> ~/.ssh/authorized_keys")
    run("chmod 600 ~/.ssh/authorized_keys")
    run("rm -f ~/id_rsa.pub")

def user_add(user, group=''):
    "Create new users on the remote server"
    password = getpass('*** enter a password for user %s:' % user)

    if group:
        sudo("useradd -m -s /bin/bash %s -g %s" % (user, group))
    else:
        sudo("useradd -m -s /bin/bash %s" % user)

    sudo('echo %s:%s | chpasswd' % (user, password))

def test():
    user_add(user="mary")

def update_distro():

    # Update the linux distribution.
    sudo("apt-get update")
    sudo("apt-get upgrade -y --show-upgraded")

    # Create group and users that will run the sertver.
    sudo("groupadd admin")

    # Add two default users
    user_add(user="www")
    user_add(user="admin", group="admin")

    # Install requirements.
    sudo("apt-get install -y postgresql postgresql-contrib postgresql-server-dev-all software-properties-common")
    sudo("apt-get install -y nginx fail2ban redis-server ufw python-software-properties g++ make")
    sudo("apt-get install -y build-essential ncurses-dev byacc zlib1g-dev python-dev git supervisor")
    sudo("apt-get install -y python-setuptools")

    # Install pip.
    sudo("easy_install pip")

    # Install the virtual environments.
    sudo("pip install virtualenv virtualenvwrapper")

    # Start webserver
    sudo("service nginx start")

    # Enable firewall.
    sudo("ufw allow ssh")
    sudo("ufw allow http")
    sudo("ufw enable")


    # Install the lessc compiler.
    sudo("sudo add-apt-repository ppa:chris-lea/node.js")
    sudo("apt-get update")
    sudo("apt-get install -y nodejs")
    sudo("npm install -g less")


