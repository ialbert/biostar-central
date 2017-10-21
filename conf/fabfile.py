from fabric.api import *
import sys

if tuple(sys.version_info) > (2, 8):
    print("Error: Fabric requires Python 2.7")
    sys.exit(-1)

main_path = '/export/sites/main_engine'
test_path = '/export/sites/test_engine'

test_env = ('source ~/miniconda3/envs/engine/bin/activate test_env')
main_env = ('source ~/miniconda3/envs/engine/bin/activate main_env')


def deploy_remote(path, env):
    with cd(path), prefix(env):
        run('git pull')
        run('python manage.py migrate')
        run('python manage.py collectstatic --noinput -v 0')
        restart_uwsgi()

def deploy_latest(path, env):
    with cd(path), prefix(env):
        run("rm -f export/engine.db")
        run('git pull')
        run('make reset')
        run('python manage.py migrate')
        run('python manage.py collectstatic --noinput -v 0')
        restart_uwsgi()

def deploy_test():
    deploy_latest(test_path, test_env)

def deploy_main():
    deploy_remote(main_path, main_env)

def restart_nginx():
    sudo("service nginx restart")

def restart_uwsgi():
    sudo("supervisorctl restart all")
