from fabric.api import run, cd
from fabric.context_managers import prefix
import socket

INSTALL_DIR = "~/sites"
INSTALL_NAME = "biostar-central"
BIOSTAR_CLONE_URL = "https://github.com/ialbert/biostar-central.git"
BRANCH = "biostar2"

BIOSTAR_HOME = "%s/%s" % (INSTALL_DIR, INSTALL_NAME)

ENV_NAME = "biostar"
LOAD_ENV = "source /usr/local/bin/virtualenvwrapper.sh"
WORKON = "%s && workon %s" % (LOAD_ENV, ENV_NAME)


def init_site():
    # Create directories.
    run('mkdir -p %s' % INSTALL_DIR)

    # Clone from repository.
    run("git clone %s %s" % (BIOSTAR_CLONE_URL, BIOSTAR_HOME))

    with cd(BIOSTAR_HOME):
        run("git checkout %s" % BRANCH)

    with prefix(LOAD_ENV):
        run("mkvirtualenv %s" % ENV_NAME)

    with prefix(WORKON):
        run("pip install -r conf/requirements/base.txt")
        run("pip install -r conf/requirements/backends.txt")
        run("pip install -r conf/requirements/celery.txt")

def init_python():
    with prefix(WORKON):
        with cd(BIOSTAR_HOME):
            run("which python")

