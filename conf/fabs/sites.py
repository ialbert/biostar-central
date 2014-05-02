from fabric.api import *

BIOSTAR_HOME = "/home/www/sites/biostar-central"
VIRTUALENV_WRAPPER = "source /usr/local/bin/virtualenvwrapper.sh"
CLONE_URL = "https://github.com/ialbert/biostar-central.git"
BRANCH = "master"

def setenv():
    # The python environment that the system needs.
    env.biostar_home = BIOSTAR_HOME
    env.wrapper = VIRTUALENV_WRAPPER
    env.biostar_clone = CLONE_URL
    env.biostar_branch = BRANCH

    # The is the main environment that will be applied to each command.
    # This is the prefix invoked when opertating on the deployed site.
    env.biostar_live = "%(biostar_home)s/live" % env
    env.biostar_env = "%(biostar_live)s/deploy.env" % env
    env.workon = "source /usr/local/bin/virtualenvwrapper.sh && workon biostar && cd %(biostar_home)s && source %(biostar_env)s" % env

def usegalaxy(user="www"):
    "Sets the environment for the biostar galaxy"
    setenv()
    env.hosts.append('biostar.usegalaxy.org')
    env.user = user

def metastars(user="www"):
    "Sets the environment for the biostar galaxy"
    setenv()
    env.hosts.append('metastars.org')
    env.user = user

def main_biostars(user="www"):
    "Sets the environment for the biostar galaxy"
    setenv()
    env.hosts.append('biostars.org')
    env.user = user

def test_site(user='www'):
    setenv()
    env.hosts.append('test.biostars.org')
    env.user = user

def hostname():
    run("hostname")
