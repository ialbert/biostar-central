from fabric.context_managers import prefix
from fabric.api import *
from fabric.contrib.files import exists
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
    env.biostar_live = "%(biostar_home)s/live/" % env
    env.workon = "source /usr/local/bin/virtualenvwrapper.sh && workon biostar && cd %(biostar_home)s && source %(biostar_env)s" % env

def server_restart():
    sudo("service nginx restart")
    sudo("supervisorctl restart biostar")

def server_config():

    with prefix(env.workon):

        # Logging into this directory.
        if not exists("live/logs"):
            run("mkdir -p live/logs")

        # Copy and link the various config files but don't overwrite
        if not exists("live/deploy.env"):
            put("conf/defaults.env", "%(biostar_live)s/deploy.env" % env)

        if not exists("live/deploy.py"):
            put("biostar/settings/deploy.py", env.biostar_live)

        if not exists("live/biostar.nginx.conf"):
            put("conf/server/biostar.nginx.conf",  env.biostar_live)
            sudo("ln -fs %(biostar_live)s/biostar.nginx.conf /etc/nginx/sites-enabled/" % env)
            sudo("service nginx restart")

        if not exists("live/biostar.supervisor.conf"):
            put("conf/server/biostar.supervisor.conf", env.biostar_live)
            sudo("ln -fs %(biostar_live)s/biostar.supervisor.conf /etc/supervisor/conf.d/" % env)
            sudo("supervisorctl restart biostar")

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

def index_biostar():
    with prefix(env.workon):
        run("./biostar.sh index")

def update_biostar():
    # Clone from repository.

    with prefix(env.workon):
        run("git pull")
        run("./biostar.sh test")
        run("python manage.py collectstatic --noinput")
