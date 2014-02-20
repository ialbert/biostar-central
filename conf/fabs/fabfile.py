from fabric.context_managers import prefix
from fabric.api import *
from getpass import getpass
from sites import *


def setenv():
    # The python environment that the system needs.
    env.biostar_home = "~/sites/biostar-central"
    env.wrapper = "source /usr/local/bin/virtualenvwrapper.sh"
    env.biostar_clone = "https://github.com/ialbert/biostar-central.git"
    env.biostar_branch = "biostar2"
    env.biostar_env = "conf/defaults.env"

    # This is the prefix invoked when opertating on the deployed site.
    env.workon = "source /usr/local/bin/virtualenvwrapper.sh && workon biostar && cd %(biostar_home)s && source %(biostar_env)s" % env

def init_config():
    
    with prefix(env.workon):
        run("cp -i conf/server/nginx.conf data/")

def init_biostar():
    # Create directories.
    run('mkdir -p %(biostar_home)s' % env)

    # Clone from repository.
    run("git clone %(biostar_clone)s %(biostar_home)s" % env)

    with cd(env.biostar_home):
        run("git checkout %(biostar_branch)s" % env)

    with prefix(env.wrapper):
        run("mkvirtualenv biostar")

    with prefix(env.workon):
        run("pip install -r %(biostar_home)s/conf/requirements/base.txt" % env)

def update_biostar():
    # Clone from repository.

    with prefix(env.workon):
        run("git pull")
        run("python manage.py collectstatic --noinput")
        sudo("sudo supervisorctl restart biostar")
        #sudo("sudo supervisorctl restart celery")
