from fabric.context_managers import prefix
from fabric.api import *
from fabric.contrib.files import exists
from getpass import getpass
from sites import *

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


def create_biostar():
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
    sudo("service nginx restart")
    sudo("supervisorctl restart biostar")
    sudo("supervisorctl restart worker beat")

def init_biostar():
    with prefix(env.workon):
        run("./biostar.sh init")
        run("./biostar.sh index")

def test():
    with prefix(env.workon):
        run("./biostar.sh test")

def migrate():
    # Clone from repository.
    with prefix(env.workon):
        run("git pull")
        run("python manage.py migrate")

def pull():
    # Perform a pull.
    with prefix(env.workon):
        run("git pull")
        #run("./biostar.sh test")
        run("python manage.py collectstatic --noinput")


