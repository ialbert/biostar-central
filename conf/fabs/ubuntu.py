"""
Command run when intializing an Ubuntu based linux distro from scratch
"""

from fabric.api import run, cd
from fabric.context_managers import prefix
from fabric.api import *
from getpass import getpass
from sites import *


def postgres_setup():
    # sudo su - postgres
    # createuser www
    pass


def add_ssh_key():
    "Appends the current SSH pub key to the remote authorized keys"
    put("~/.ssh/id_rsa.pub", "~/")
    run("mkdir -p .ssh")
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
    # Set the hostname
    host = prompt('enter the hostname')
    ip = prompt('enter the ip number')

    sudo('echo  %s > /etc/hostname' % host)
    sudo('hostname -F /etc/hostname')
    sudo('echo "127.0.0.1     localhost.localdomain    localhost" > /etc/hosts')
    sudo('echo "%s    %s    %s" >> /etc/hosts' % ip, host, host)

    # Update the linux distribution.
    sudo("apt-get update")
    sudo("apt-get upgrade -y --show-upgraded")

    # Create group and users that will run the sertver.
    sudo("groupadd admin")

    # Install requirements.
    sudo("apt-get install -y postgresql postgresql-contrib postgresql-server-dev-all software-properties-common")
    sudo("apt-get install -y nginx fail2ban redis-server ufw python-software-properties g++ make openjdk-7-jdk")
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

    # Installing elastic search
    # wget https://download.elasticsearch.org/elasticsearch/elasticsearch/elasticsearch-1.0.0.deb
    # sudo dpkg -i elasticsearch-1.0.1.deb
    # sudo update-rc.d elasticsearch defaults 95 10
    # sudo /etc/init.d/elasticsearch start




def install_nodejs():

    # reconfigure timezone
    # dpkg-reconfigure tzdata

    # Add two default users
    user_add(user="www", group="admin")

    # Install the lessc compiler.
    sudo("ufw enable")
    sudo("sudo add-apt-repository ppa:chris-lea/node.js")
    sudo("apt-get install -y nodejs")
    sudo("npm install -g less")


