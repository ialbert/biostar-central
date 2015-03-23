"""
Fabric command manager.
"""
from fabric.context_managers import prefix
from fabric.api import *
from fabric.contrib.files import exists
from getpass import getpass
from fabric.api import *

BIOSTAR_HOME = "/export/sites/{}.biostars.org"
VIRTUALENV_WRAPPER = "source /usr/local/bin/virtualenvwrapper.sh"
CLONE_URL = "https://github.com/ialbert/biostar-central.git"
BRANCH = "master"

def setenv(domain):
    # The python environment that the system needs.
    env.biostar_home = BIOSTAR_HOME.format(domain)
    env.wrapper = VIRTUALENV_WRAPPER
    env.biostar_clone = CLONE_URL
    env.biostar_branch = BRANCH

    # The is the main environment that will be applied to each command.
    # This is the prefix invoked when opertating on the deployed site.
    env.biostar_live = "%(biostar_home)s/live" % env
    env.biostar_env = "%(biostar_live)s/deploy.env" % env
    env.workon = "source /usr/local/bin/virtualenvwrapper.sh && workon b3 && cd %(biostar_home)s && source %(biostar_env)s" % env


def test_site(user='www'):
    setenv(domain="test")
    env.hosts.append('www.test.biostars.org')
    env.user = user

def hostname():
    """
    Runs a test command
    '"""
    run("hostname")

def pull():
    """
    Perform a pull followe by a biostar init on the remote site.
    """
    with prefix(env.workon):
        run("git pull")
        run("python manage.py migrate")
        run("python manage.py collectstatic --noinput")

