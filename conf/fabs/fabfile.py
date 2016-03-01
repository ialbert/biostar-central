from fabric.context_managers import prefix
from fabric.api import *
from fabric.contrib.files import exists
from getpass import getpass

BIOSTAR_HOME = "/home/www/sites/biostar-central"
VIRTUALENV_WRAPPER = "source /usr/local/bin/virtualenvwrapper.sh"
CLONE_URL = "https://github.com/ialbert/biostar-central.git"
BRANCH = "production"

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

def biostars(user="www"):
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

def copy_config():
    # Create a default environment.
    if not exists(env.biostar_env):
        put("conf/defaults.env", env.biostar_env)

        # Logging into this directory.
    if not exists("%(biostar_live)s/logs" % env):
        run("mkdir -p %(biostar_live)s/logs" % env)

    # Customize the deployment environment.
    with prefix(env.workon):

        # Copy over all scripts
        scripts = [
            "celery.beat.sh",
            "celery.worker.sh",
            "gunicorn.start.sh",
        ]
        for name in scripts:
            if not exists("live/%s" % name):
                put("conf/server/%s" % name, env.biostar_live)
                run("chmod +x %(biostar_live)s/*.sh" % env)

        if not exists("live/deploy.py"):
            put("biostar/settings/deploy.py", env.biostar_live)

        if not exists("live/biostar.nginx.conf"):
            put("conf/server/biostar.nginx.conf", env.biostar_live)
            sudo("ln -fs %(biostar_live)s/biostar.nginx.conf /etc/nginx/sites-enabled/" % env)

        if not exists("live/biostar.supervisor.conf"):
            put("conf/server/biostar.supervisor.conf", env.biostar_live)
            sudo("ln -fs %(biostar_live)s/biostar.supervisor.conf /etc/supervisor/conf.d/" % env)


def create_directories():
    "Create the biostar directories"

    # Create directories.
    run('mkdir -p %(biostar_home)s' % env)

    # Clone from repository.
    run("git clone %(biostar_clone)s %(biostar_home)s" % env)

    with cd(env.biostar_home):
        run("git checkout %(biostar_branch)s" % env)

        with prefix(env.wrapper):
            run("mkvirtualenv biostar")
            run("workon biostar && pip install -r %(biostar_home)s/conf/requirements/all.txt" % env)


def restart():
    """A full restart"""
    sudo("service nginx restart")
    sudo("supervisorctl restart biostar")
    sudo("supervisorctl restart worker beat")


def init():
    """Initialize the website"""
    with prefix(env.workon):
        run("./biostar.sh init")


def index():
    """Index the content"""
    with prefix(env.workon):
        run("./biostar.sh index")


def test():
    """Run test"""
    with prefix(env.workon):
        run("./biostar.sh test")


def migrate():
    """Migrate from main"""
    with prefix(env.workon):
        run("python manage.py migrate")


def pull():
    """Perform a pull on remote"""
    with prefix(env.workon):
        run("git pull")

def pip():
     """Reinstall requirments"""
     with prefix(env.workon):
        run("pip install -r conf/requirements/base.txt")


def deploy():
    """A full deploy"""
    pull()
    pip()
    init()
    restart()
    