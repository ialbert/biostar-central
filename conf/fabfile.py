from fabric.api import *
import sys

if tuple(sys.version_info) > (2, 8):
    print("Error: Fabric requires Python 2.7")
    sys.exit(-1)

main_path = '/export/sites/main_engine'
test_path = '/export/sites/test_engine'

env_prefix = ('source ~/miniconda3/envs/engine/bin/activate engine')


def deploy_remote(path):
    with cd(path), prefix(env_prefix):

        run('git pull')
        run('python manage.py migrate')
        run('python manage.py collectstatic --noinput -v 0')


def deploy_latest(path):
    with cd(path), prefix(env_prefix):

        run("rm -f export/engine.db")
        run('git pull https://github.com/Natay/biostar-engine.git')
        run('python manage.py migrate')
        run('python manage.py collectstatic --noinput -v 0')


def deploy_test():
    deploy_latest(test_path)

def deploy_main():
    deploy_remote(main_path)

def restart_nginx():
    sudo("service nginx restart")

def restart_uwsgi():
    sudo("supervisorctl restart all")